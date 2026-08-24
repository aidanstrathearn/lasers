use crate::LaserApp;
use crate::plotter::Plotter;
use laser_solver::error::SolverError;

impl LaserApp {
    pub fn kappa_plot(&self) -> Result<Plotter, SolverError> {
        let z = self.grid_points.grid(self.fibre_params.length);
        let kappas = self.grating.grid(self.grid_points.0);

        // kappas.len() = z.len() - 1, so zip ignores last element of z
        let points = z
            .into_iter()
            .zip(kappas.into_iter())
            .map(|(z, k)| [z, k])
            .collect();

        let mut plot = Plotter::new();
        plot.add_points(points).label("Kappa");
        plot.xlabel("z");
        plot.ylabel("Kappa");
        Ok(plot)
    }
}
