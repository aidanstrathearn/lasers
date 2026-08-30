pub mod picard;
mod shooting;

use crate::dopant::{DopantModel, TwoLevelDopant};
use crate::error::SolverError;
use crate::grating::{GratingModel, PiShift, sample_grating};
use crate::lase::{
    FieldProfile, Pump, PumpScan, ResolvedFibre, UniformGrid,
    find_threshold_and_slope as scan_for_threshold, pump_scan as scan_pump_totals,
};
use crate::maths::picard::PicardConfig;
use crate::maths::rootfind::RootFindConfig;
use crate::maths::utils::IterationConfig;

#[derive(Copy, Clone)]
pub struct DfbSolveConfig {
    pub root_find: RootFindConfig,
    pub picard: PicardConfig,
}

pub struct DfbLaser<
    'a,
    D: DopantModel = TwoLevelDopant,
    G: GratingModel = PiShift,
> {
    pub fibre: ResolvedFibre<'a, D, G>,
    grid: UniformGrid,
    kappas: Box<[f64]>,
}

impl<'a, D: DopantModel, G: GratingModel> DfbLaser<'a, D, G> {
    pub fn new(fibre: ResolvedFibre<'a, D, G>, steps: usize) -> Self {
        let grid = UniformGrid::new(fibre.length(), steps);
        let kappas = sample_grating(fibre.grating(), grid.steps()).into_boxed_slice();
        DfbLaser {
            fibre,
            grid,
            kappas,
        }
    }

    pub fn grid(&self) -> UniformGrid {
        self.grid
    }

    pub fn kappas(&self) -> &[f64] {
        &self.kappas
    }

    pub fn solve(
        &self,
        pump: Pump,
        solve_config: DfbSolveConfig,
        full_profile: bool,
    ) -> Result<FieldProfile, SolverError> {
        let use_picard = pump.backward_amplitude() > 0.0;
        if use_picard {
            self.solve_picard(pump, solve_config, full_profile)
        } else {
            self.solve_shooting(pump, solve_config, full_profile)
        }
    }

    pub fn find_threshold_and_slope(
        &self,
        pump_start: Pump,
        pump_step: f64,
        iteration: IterationConfig,
        solve_config: DfbSolveConfig,
    ) -> Result<(f64, f64, f64), SolverError> {
        pump_start.amplitudes();
        let use_picard = pump_start.balance != 1.0;
        if use_picard {
            let mut solver = self.initial_picard_solver(pump_start);
            let f = |total| {
                self.output_power_picard(
                    Pump {
                        total,
                        ..pump_start
                    },
                    solve_config,
                    &mut solver,
                )
            };
            scan_for_threshold(pump_start.total, pump_step, iteration, f)
        } else {
            let f = |total| {
                self.output_power_shooting(
                    Pump {
                        total,
                        ..pump_start
                    },
                    solve_config,
                )
            };
            scan_for_threshold(pump_start.total, pump_step, iteration, f)
        }
    }

    pub fn pump_scan(
        &self,
        pump_totals: &[f64],
        balance: f64,
        solve_config: DfbSolveConfig,
    ) -> Result<PumpScan, SolverError> {
        let Some(&pump_start) = pump_totals.first() else {
            return Ok(Vec::new());
        };

        let use_picard = balance != 1.0;
        if use_picard {
            let mut solver = self.initial_picard_solver(Pump {
                total: pump_start,
                balance,
            });
            scan_pump_totals(pump_totals, |total| {
                self.output_power_picard(Pump { total, balance }, solve_config, &mut solver)
            })
        } else {
            scan_pump_totals(pump_totals, |total| {
                self.output_power_shooting(Pump { total, balance }, solve_config)
            })
        }
    }
}
