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
        set_boundary: SetBoundary,
        step: Step,
        error: Error,
    ) -> Result<&[T], PicardError>
    where
        SetBoundary: FnMut(&[T]) -> T,
        Step: FnMut(&T, &T, usize) -> T,
        Error: FnMut(&[T], &[T]) -> f64,
    {
        self.solve_relaxed(max_iterations, set_boundary, step, error, |_, _| {})
    }

    pub fn solve_relaxed<SetBoundary, Step, Error, Relax>(
        &mut self,
        max_iterations: usize,
        mut set_boundary: SetBoundary,
        mut step: Step,
        mut error: Error,
        mut relax: Relax,
    ) -> Result<&[T], PicardError>
    where
        SetBoundary: FnMut(&[T]) -> T,
        Step: FnMut(&T, &T, usize) -> T,
        Error: FnMut(&[T], &[T]) -> f64,
        Relax: FnMut(&mut T, &T),
    {
        for _ in 0..max_iterations {
            self.new[0] = set_boundary(&self.current);
            for i in 0..self.new.len() - 1 {
                self.new[i + 1] = step(&self.new[i], &self.current[i], i);
            }

            if error(&self.new, &self.current) <= 1.0 {
                std::mem::swap(&mut self.current, &mut self.new);
                return Ok(&self.current);
            }

            for (new, current) in self.new.iter_mut().zip(&self.current) {
                relax(new, current);
            }

            std::mem::swap(&mut self.current, &mut self.new);
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

    #[test]
    fn relaxed_solve_blends_raw_profile_with_current_profile_in_place() {
        let mut solver = PicardSolver::from_initial(vec![0.0, 0.0]);
        let result = solver.solve_relaxed(
            2,
            |_| 4.0,
            |new_previous, _, _| new_previous + 4.0,
            |_, _| 2.0,
            |new, current| *new = 0.5 * (*new + *current),
        );

        assert!(matches!(result, Err(PicardError::DidNotConverge)));
        assert_eq!(solver.profile(), [3.0, 6.0]);
    }

    #[test]
    fn relaxed_solve_checks_raw_profile_before_relaxing() {
        let mut solver = PicardSolver::from_initial(vec![0.0_f64]);
        let mut relaxation_calls = 0;

        let result = solver.solve_relaxed(
            1,
            |_| 1.0,
            |_, _, _| unreachable!(),
            |new, current| (new[0] - current[0]).abs(),
            |_, _| relaxation_calls += 1,
        );

        assert!(result.is_ok());
        assert_eq!(relaxation_calls, 0);
        assert_eq!(solver.profile(), [1.0]);
    }
}
