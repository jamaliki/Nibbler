# Security policy

## Supported versions

Until the first stable release, security fixes target the latest commit on `main` and
the newest published `0.x` release, if one exists. Older alpha builds are not supported.

## Reporting a vulnerability

Please report suspected vulnerabilities through a
[private GitHub security advisory](https://github.com/jamaliki/Nibbler/security/advisories/new).
Do not include exploit details in a public issue.

Include the affected Nibbler version or commit, platform, input format, minimal
reproduction, observed impact, and whether untrusted input is required. Parser crashes,
unbounded resource use, validation bypasses, unsafe-code concerns, and inconsistent
resource-limit enforcement are all appropriate reports.

The project will acknowledge the report, investigate it privately, and coordinate a
fix and disclosure according to severity. Please allow a reasonable remediation period
before public disclosure.

Nibbler parses untrusted data without runtime network access, but no software should be
treated as a security boundary unless its deployment has been assessed for that use.
