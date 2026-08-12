use std::collections::BTreeMap;
use std::io;
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::sync::{Arc, Mutex, MutexGuard};
use std::thread::{self, JoinHandle};

use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyModule};

use crate::cif::ProjectionPlan;
use crate::python::{
    ErrorFields, InputSource, PredicateSpec, ReadFailure, ReadOutput, build_optional_plan,
    error_fields, error_to_python, execute_read, output_to_python,
};

type Job = (usize, InputSource);
type Completed = (usize, Result<ReadOutput, ReadFailure>);
type PythonScanItem = (usize, Option<Py<PyAny>>, Option<ErrorFields>);

#[pyclass(name = "_NativeScan")]
struct NativeScan {
    job_sender: Option<SyncSender<Job>>,
    completed_receiver: Arc<Mutex<Receiver<Completed>>>,
    workers: Vec<JoinHandle<()>>,
    completed_by_index: BTreeMap<usize, Result<ReadOutput, ReadFailure>>,
    submitted_count: usize,
    yielded_count: usize,
    schema: Option<String>,
}

#[pymethods]
impl NativeScan {
    #[new]
    #[pyo3(signature = (workers, category=None, columns=None, predicates=Vec::new(), schema=None))]
    fn new(
        workers: usize,
        category: Option<String>,
        columns: Option<Vec<String>>,
        predicates: Vec<PredicateSpec>,
        schema: Option<String>,
    ) -> PyResult<Self> {
        if workers == 0 {
            return Err(PyValueError::new_err("workers must be greater than zero"));
        }
        let plan = build_optional_plan(category, columns, predicates, schema.clone())
            .map_err(error_to_python)?;
        Self::spawn(workers, plan, schema)
            .map_err(|error| PyRuntimeError::new_err(format!("cannot start scan worker: {error}")))
    }

    fn submit_file(&mut self, py: Python<'_>, file: String) -> PyResult<()> {
        self.submit(py, InputSource::File(file))
    }

    fn submit_bytes(
        &mut self,
        py: Python<'_>,
        source_name: String,
        bytes: Vec<u8>,
    ) -> PyResult<()> {
        self.submit(py, InputSource::Bytes { source_name, bytes })
    }

    fn next<'py>(&mut self, py: Python<'py>) -> PyResult<Option<PythonScanItem>> {
        if self.yielded_count == self.submitted_count {
            return Ok(None);
        }
        let index = self.yielded_count;
        while !self.completed_by_index.contains_key(&index) {
            let receiver = Arc::clone(&self.completed_receiver);
            let completed = py
                .detach(move || lock_receiver(&receiver).recv())
                .map_err(|_| PyRuntimeError::new_err("native scan workers stopped early"))?;
            self.completed_by_index.insert(completed.0, completed.1);
        }
        let Some(result) = self.completed_by_index.remove(&index) else {
            return Err(PyRuntimeError::new_err(
                "ordered scan result disappeared from its reorder buffer",
            ));
        };
        self.yielded_count += 1;
        match result {
            Ok(output) => Ok(Some((
                index,
                Some(output_to_python(py, output, true, self.schema.clone())?),
                None,
            ))),
            Err(error) => Ok(Some((index, None, Some(error_fields(&error))))),
        }
    }

    fn close(&mut self, py: Python<'_>) {
        self.job_sender.take();
        let workers = std::mem::take(&mut self.workers);
        py.detach(move || {
            for worker in workers {
                let _ = worker.join();
            }
        });
    }

    fn cancel(&mut self) {
        self.job_sender.take();
        self.workers.clear();
    }
}

impl NativeScan {
    fn spawn(
        worker_count: usize,
        plan: Option<ProjectionPlan>,
        schema: Option<String>,
    ) -> io::Result<Self> {
        let (job_sender, job_receiver) = mpsc::sync_channel(worker_count);
        let job_receiver = Arc::new(Mutex::new(job_receiver));
        let (completed_sender, completed_receiver) = mpsc::channel();
        let mut workers = Vec::with_capacity(worker_count);
        for worker_index in 0..worker_count {
            workers.push(spawn_worker(
                worker_index,
                Arc::clone(&job_receiver),
                completed_sender.clone(),
                plan.clone(),
            )?);
        }
        drop(completed_sender);
        Ok(Self {
            job_sender: Some(job_sender),
            completed_receiver: Arc::new(Mutex::new(completed_receiver)),
            workers,
            completed_by_index: BTreeMap::new(),
            submitted_count: 0,
            yielded_count: 0,
            schema,
        })
    }

    fn submit(&mut self, py: Python<'_>, input: InputSource) -> PyResult<()> {
        let Some(sender) = &self.job_sender else {
            return Err(PyRuntimeError::new_err("cannot submit to a closed scan"));
        };
        let sender = sender.clone();
        let index = self.submitted_count;
        py.detach(move || sender.send((index, input)))
            .map_err(|_| PyRuntimeError::new_err("native scan workers stopped early"))?;
        self.submitted_count += 1;
        Ok(())
    }
}

impl Drop for NativeScan {
    fn drop(&mut self) {
        self.job_sender.take();
        // Workers hold no Python objects and exit after already-admitted jobs. Detaching
        // here ensures Python garbage collection never blocks while holding the GIL.
        self.workers.clear();
    }
}

fn spawn_worker(
    worker_index: usize,
    job_receiver: Arc<Mutex<Receiver<Job>>>,
    completed_sender: mpsc::Sender<Completed>,
    plan: Option<ProjectionPlan>,
) -> io::Result<JoinHandle<()>> {
    thread::Builder::new()
        .name(format!("nibbler-scan-{worker_index}"))
        .spawn(move || {
            loop {
                let job = lock_receiver(&job_receiver).recv();
                let Ok((index, input)) = job else {
                    break;
                };
                if completed_sender
                    .send((index, execute_read(input, plan.clone())))
                    .is_err()
                {
                    break;
                }
            }
        })
}

fn lock_receiver<T>(receiver: &Arc<Mutex<Receiver<T>>>) -> MutexGuard<'_, Receiver<T>> {
    receiver
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<NativeScan>()?;
    Ok(())
}
