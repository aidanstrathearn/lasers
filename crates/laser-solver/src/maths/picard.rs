use std::fmt;

#[derive(Debug)]
pub enum PicardError {
    DidNotConverge,
}

impl fmt::Display for PicardError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DidNotConverge => {
                write!(formatter, "picard iteration did not converge")
            }
        }
    }
}

impl std::error::Error for PicardError {}

#[derive(Copy, Clone, Debug)]
pub struct PicardConfig {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
}

impl Default for PicardConfig {
    fn default() -> Self {
        Self {
            max_iterations: 500,
            relative_tolerance: 1e-8,
            absolute_tolerance: 1e-12,
        }
    }
}

pub struct PicardSolver<T> {
    current: Vec<T>,
    new: Vec<T>,
}

impl<T: Clone> PicardSolver<T> {
    pub fn filled(len: usize, initial: T) -> Self {
        Self::from_initial(vec![initial; len])
    }

    pub fn from_initial(initial: Vec<T>) -> Self {
        assert!(
            !initial.is_empty(),
            "Picard initial profile must not be empty"
        );
        Self {
            new: initial.clone(),
            current: initial,
        }
    }
}

impl<T> PicardSolver<T> {
    pub fn profile(&self) -> &[T] {
        &self.current
    }

    pub fn solve<SetBoundary, Step, Error>(
        &mut self,
        max_iterations: usize,
        mut set_boundary: SetBoundary,
        mut step: Step,
        mut error: Error,
    ) -> Result<&[T], PicardError>
    where
        SetBoundary: FnMut(&[T]) -> T,
        Step: FnMut(&T, &T, usize) -> T,
        Error: FnMut(&[T], &[T]) -> f64,
    {
        for _ in 0..max_iterations {
            self.new[0] = set_boundary(&self.current);
            for i in 0..self.new.len() - 1 {
                self.new[i + 1] = step(&self.new[i], &self.current[i], i);
            }

            std::mem::swap(&mut self.current, &mut self.new);
            if error(&self.current, &self.new) <= 1.0 {
                return Ok(&self.current);
            }
        }
        Err(PicardError::DidNotConverge)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_a_constant_initial_profile() {
        let solver = PicardSolver::filled(3, "initial".to_owned());

        assert_eq!(solver.profile(), ["initial", "initial", "initial"]);
    }

    #[test]
    fn step_receives_new_previous_old_current_and_interval_index() {
        let initial = vec![1.0, 2.0, 3.0];
        let mut solver = PicardSolver::from_initial(initial);
        let mut steps = Vec::new();

        solver
            .solve(
                1,
                |current| {
                    assert_eq!(current, [1.0, 2.0, 3.0]);
                    10.0
                },
                |new_previous, old_current, i| {
                    steps.push((*new_previous, *old_current, i));
                    new_previous + 1.0
                },
                |current, previous| {
                    assert_eq!(current, [10.0, 11.0, 12.0]);
                    assert_eq!(previous, [1.0, 2.0, 3.0]);
                    0.0
                },
            )
            .expect("error callback should accept the first Picard iteration");

        assert_eq!(steps, vec![(10.0, 1.0, 0), (11.0, 2.0, 1)]);
        assert_eq!(solver.profile(), [10.0, 11.0, 12.0]);
    }

    #[test]
    #[should_panic(expected = "Picard initial profile must not be empty")]
    fn rejects_an_empty_initial_profile() {
        let _solver = PicardSolver::<String>::from_initial(Vec::new());
    }

    #[test]
    fn supports_non_copy_states() {
        let mut solver = PicardSolver::from_initial(vec!["old 0".to_owned(), "old 1".to_owned()]);

        solver
            .solve(
                1,
                |_| "new 0".to_owned(),
                |new_previous, old_current, _| format!("{new_previous}; {old_current}"),
                |_, _| 0.0,
            )
            .unwrap();

        assert_eq!(solver.profile(), ["new 0", "new 0; old 0"]);
    }
}
