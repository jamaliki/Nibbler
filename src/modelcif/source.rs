//! Construction of typed ModelCIF metadata over one PDBx coordinate model.

use crate::cif::{CategoryView, CifDocument};
use crate::pdbx::{ComponentRegistry, SemanticError, build_model_with_registry as build_pdbx};

use super::aggregate::ModelCifModel;
use super::{prediction, provenance, qa, template};

/// Build ModelCIF semantics without a caller-selected component registry.
///
/// # Errors
///
/// Returns an error when the coordinate graph or a present ModelCIF row is structurally invalid.
pub fn build_model(document: &CifDocument) -> Result<ModelCifModel, SemanticError> {
    build_model_with_registry(document, None)
}

/// Build ModelCIF semantics with an optional immutable local CCD registry.
///
/// # Errors
///
/// Returns an error when the coordinate graph or a present ModelCIF row is structurally invalid.
pub fn build_model_with_registry(
    document: &CifDocument,
    registry: Option<&ComponentRegistry>,
) -> Result<ModelCifModel, SemanticError> {
    let coordinates = build_pdbx(document, registry)?;
    let categories = CategoryView::new(document.blocks()[0].entries());

    Ok(ModelCifModel {
        coordinates,
        audit_conform: prediction::audit_conform(categories)?,
        data: prediction::data(categories)?,
        data_groups: prediction::data_groups(categories)?,
        targets: prediction::targets(categories)?,
        target_instances: prediction::target_instances(categories)?,
        models: prediction::models(categories)?,
        model_groups: prediction::model_groups(categories)?,
        model_group_links: prediction::model_group_links(categories)?,
        representatives: prediction::representatives(categories)?,
        software: provenance::software(categories)?,
        software_groups: provenance::software_groups(categories)?,
        protocol_steps: provenance::protocol_steps(categories)?,
        templates: template::templates(categories)?,
        template_segments: template::template_segments(categories)?,
        template_mappings: template::template_mappings(categories)?,
        alignments: template::alignments(categories)?,
        alignment_details: template::alignment_details(categories)?,
        alignment_sequences: template::alignment_sequences(categories)?,
        qa_metrics: qa::qa_metrics(categories)?,
        qa_values: qa::qa_values(categories)?,
        associated_files: provenance::associated_files(categories)?,
        archive_members: provenance::archive_members(categories)?,
    })
}

#[cfg(test)]
#[allow(clippy::expect_used)]
mod tests {
    use crate::cif::parse;

    use super::super::qa::QaValue;
    use super::*;

    const PREDICTION: &[u8] =
        include_bytes!("../../tests/fixtures/modelcif/prediction_with_qa.cif");

    #[test]
    fn retains_prediction_provenance_template_alignment_and_qa_records() {
        let source = String::from_utf8(PREDICTION.to_vec())
            .expect("fixture is UTF-8")
            .replace(
                "_software.description\n1 Rosetta 'model building' 2.1.0 'deterministic test predictor'",
                "_software.description\n_software.location\n1 Rosetta 'model building' 2.1.0 'deterministic test predictor' local://rosetta",
            )
            .replace("1 1 1 A . . . .", "1 1 1 A RMSD 1.25 50.0 aligned");
        let document = parse(source.as_bytes()).expect("fixture parses");
        let model = build_model(&document).expect("fixture builds");

        assert_eq!(
            model.audit_conform,
            [
                ("mmcif_pdbx.dic".to_owned(), "5.416".to_owned()),
                ("mmcif_ma.dic".to_owned(), "1.4.9".to_owned()),
            ]
        );
        let data = &model.data[0];
        assert_eq!(
            (data.id, data.name.as_str(), data.content_type.as_str()),
            (1, "designed target sequence", "target")
        );
        let data_group = &model.data_groups[0];
        assert_eq!(
            (
                data_group.ordinal_id,
                data_group.group_id,
                data_group.data_id
            ),
            (1, 1, 1)
        );
        let target = &model.targets[0];
        assert_eq!(
            (
                target.entity_id.as_str(),
                target.data_id,
                target.origin.as_str()
            ),
            ("1", 1, "designed")
        );
        let instance = &model.target_instances[0];
        assert_eq!(
            (
                instance.asym_id.as_str(),
                instance.entity_id.as_str(),
                instance.details.as_deref()
            ),
            ("A", "1", Some("designed peptide target"))
        );

        let prediction = &model.models[0];
        assert_eq!(
            (
                prediction.id,
                prediction.name.as_deref(),
                prediction.assembly_id,
                prediction.model_type.as_str(),
                prediction.type_details.as_deref(),
                prediction.data_id,
            ),
            (
                1,
                Some("ranked model 1"),
                Some(1),
                "Ab initio model",
                None,
                2
            )
        );
        let group = &model.model_groups[0];
        assert_eq!(
            (group.id, group.name.as_deref(), group.details.as_deref()),
            (
                1,
                Some("primary prediction"),
                Some("single-model prediction group")
            )
        );
        let link = &model.model_group_links[0];
        assert_eq!((link.model_id, link.group_id), (1, 1));
        let representative = &model.representatives[0];
        assert_eq!(
            (
                representative.id,
                representative.group_id,
                representative.model_id,
                representative.selection_criteria.as_str()
            ),
            (1, 1, 1, "best scoring model")
        );

        let software = &model.software[0];
        assert_eq!(
            (
                software.id,
                software.name.as_str(),
                software.classification.as_str(),
                software.version.as_deref(),
                software.description.as_deref(),
                software.location.as_deref(),
            ),
            (
                1,
                "Rosetta",
                "model building",
                Some("2.1.0"),
                Some("deterministic test predictor"),
                Some("local://rosetta"),
            )
        );
        let software_group = &model.software_groups[0];
        assert_eq!(
            (
                software_group.ordinal_id,
                software_group.group_id,
                software_group.software_id,
                software_group.parameter_group_id,
            ),
            (1, 1, 1, None)
        );
        let step = &model.protocol_steps[0];
        assert_eq!(
            (
                step.ordinal_id,
                step.protocol_id,
                step.step_id,
                step.method_type.as_str(),
                step.name.as_deref(),
                step.details.as_deref(),
                step.software_group_id,
                step.input_data_group_id,
                step.output_data_group_id,
            ),
            (
                1,
                1,
                1,
                "modeling",
                Some("prediction"),
                Some("coordinates generated from supplied target sequence"),
                Some(1),
                Some(1),
                Some(2),
            )
        );

        let template = &model.templates[0];
        assert_eq!(
            (
                template.ordinal_id,
                template.id,
                template.origin.as_str(),
                template.entity_type.as_str(),
                template.data_id,
                template.target_asym_id.as_str(),
                template.auth_asym_id.as_str(),
                template.label_asym_id.as_deref(),
                template.label_entity_id.as_deref(),
                template.model_number,
                template.transform_id,
                template.name.as_deref(),
            ),
            (
                1,
                1,
                "customized",
                "polymer",
                4,
                "A",
                "T",
                Some("T"),
                Some("1"),
                1,
                1,
                Some("custom test template"),
            )
        );
        let segment = &model.template_segments[0];
        assert_eq!(
            (
                segment.id,
                segment.template_id,
                segment.sequence_begin,
                segment.sequence_end
            ),
            (1, 1, Some(1), Some(2))
        );
        let mapping = &model.template_mappings[0];
        assert_eq!(
            (
                mapping.id,
                mapping.template_segment_id,
                mapping.target_asym_id.as_str(),
                mapping.target_begin,
                mapping.target_end,
            ),
            (1, 1, "A", Some(1), Some(2))
        );
        let alignment = &model.alignments[0];
        assert_eq!(
            (
                alignment.id,
                alignment.data_id,
                alignment.software_group_id,
                alignment.length,
                alignment.alignment_type.as_deref(),
                alignment.mode.as_deref(),
            ),
            (
                1,
                5,
                Some(1),
                Some(2),
                Some("target-template pairwise alignment"),
                Some("global"),
            )
        );
        let detail = &model.alignment_details[0];
        assert_eq!(
            (
                detail.ordinal_id,
                detail.alignment_id,
                detail.template_segment_id,
                detail.target_asym_id.as_str(),
                detail.score_type.as_deref(),
                detail.score_value,
                detail.sequence_identity,
            ),
            (1, 1, 1, "A", Some("RMSD"), Some(1.25), Some(50.0))
        );
        assert_eq!(
            model
                .alignment_sequences
                .iter()
                .map(|sequence| (
                    sequence.ordinal_id,
                    sequence.alignment_id,
                    sequence.target_template_flag.as_str(),
                    sequence.sequence.as_str(),
                ))
                .collect::<Vec<_>>(),
            [(1, 1, "1", "AM"), (2, 1, "2", "AM")]
        );

        let metric = &model.qa_metrics[0];
        assert_eq!(
            (
                metric.id,
                metric.name.as_str(),
                metric.description.as_deref(),
                metric.metric_type.as_str(),
                metric.mode.as_str(),
                metric.software_group_id,
                metric.data_id,
            ),
            (
                1,
                "pLDDT",
                Some("per-residue prediction confidence, not displacement"),
                "pLDDT",
                "local",
                Some(1),
                Some(3),
            )
        );
        assert!(matches!(
            &model.qa_values[0],
            QaValue::Global {
                ordinal_id: 1,
                model_id: 1,
                metric_id: 2,
                value: 0.82,
            }
        ));
        assert!(matches!(
            &model.qa_values[1],
            QaValue::Local {
                ordinal_id: 1,
                model_id: 1,
                site,
                metric_id: 1,
                value: 91.0,
            } if site.asym_id == "A" && site.sequence_id == 1 && site.component_id == "ALA"
        ));
        assert!(matches!(
            &model.qa_values[3],
            QaValue::Pairwise {
                ordinal_id: 1,
                model_id: 1,
                first,
                second,
                metric_id: 3,
                value: 2.5,
            } if first.asym_id == "A"
                && first.sequence_id == 1
                && first.component_id == "ALA"
                && second.asym_id == "A"
                && second.sequence_id == 2
                && second.component_id == "MSE"
        ));

        let associated = &model.associated_files[0];
        assert_eq!(
            (
                associated.id,
                associated.entry_id.as_str(),
                associated.file_url.as_str(),
                associated.file_type.as_deref(),
                associated.file_format.as_deref(),
                associated.file_content.as_deref(),
                associated.details.as_deref(),
                associated.data_id,
            ),
            (
                1,
                "NIBBLER_LIGAND_ION_WATER",
                "https://example.invalid/nibbler-qa.zip",
                Some("archive"),
                Some("zip"),
                Some("archive with multiple files"),
                Some("test archive; not fetched"),
                Some(3),
            )
        );
        let member = &model.archive_members[0];
        assert_eq!(
            (
                member.id,
                member.archive_file_id,
                member.file_path.as_str(),
                member.file_format.as_deref(),
                member.file_content.as_deref(),
                member.description.as_deref(),
                member.data_id,
            ),
            (
                1,
                1,
                "scores/pae.json",
                Some("json"),
                Some("QA metrics"),
                Some("pairwise QA matrix"),
                Some(3),
            )
        );
    }
}
