"""Exception hierarchy with a defined caller response for each class.

See docs/05-software-architecture.md section 5.
"""


class VirotraceError(Exception):
    """Base class for all Virotrace errors."""


class ConfigError(VirotraceError):
    """Configuration is invalid or internally inconsistent.

    Response: abort immediately, naming the offending key and the constraint.
    """


class ContractError(VirotraceError):
    """A data contract was violated (schema, coordinates, reading frame).

    Response: abort the stage. This indicates a bug or bad reference data,
    not bad input -- continuing would produce results nothing downstream can trust.
    """


class RecordError(VirotraceError):
    """A single record is unusable.

    Response: quarantine the record, increment the counter, continue. Escalates to
    an abort if the quarantine rate exceeds ``qc.max_record_error_rate`` -- silent
    mass-quarantine is the failure mode this split would otherwise invite.
    """

    def __init__(self, message: str, *, record_id: str | None = None, rule: str | None = None):
        super().__init__(message)
        self.record_id = record_id
        self.rule = rule
