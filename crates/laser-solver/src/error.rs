use crate::picard::PicardError;
use crate::rootfind::RootFindError;
use std::fmt;

#[derive(Debug)]
pub enum SolverError {
    RootFind(RootFindError),
    Picard(PicardError),
    ThresholdNotFound,
}

impl From<RootFindError> for SolverError {
    fn from(error: RootFindError) -> Self {
        Self::RootFind(error)
    }
}

impl From<PicardError> for SolverError {
    fn from(error: PicardError) -> Self {
        Self::Picard(error)
    }
}

impl fmt::Display for SolverError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RootFind(error) => write!(formatter, "{error}"),
            Self::Picard(error) => write!(formatter, "{error}"),
            Self::ThresholdNotFound => write!(formatter, "threshold not found"),
        }
    }
}

impl std::error::Error for SolverError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::RootFind(error) => Some(error),
            Self::Picard(error) => Some(error),
            Self::ThresholdNotFound => None,
        }
    }
}
