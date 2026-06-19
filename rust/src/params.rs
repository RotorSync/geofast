//! AircraftParams — calibrated constants. Direct port of the Python dataclass.
//! All distances in feet, speeds ft/s, times seconds.

pub const ACRE_FT2: f64 = 43560.0;

#[derive(Debug, Clone, Copy)]
pub struct AircraftParams {
    pub swath_ft: f64,     // effective swath (R44 booms, 50-55 ft real)
    pub speed_fps: f64,    // 63 mph (calibrated from 49k-ac sample)
    pub t_turn_s: f64,     // avg ag turnaround (median of 15803 measured)
    pub tank_gal: f64,     // effective usable load (~65 of 80 gal usable)
    pub gpa: f64,          // default gpa; per-job rate parsed when present
    pub t_load_s: f64,     // hot load median (mean 113s, median 100s, 1853 loads)
    pub t_cycle_s: f64,    // per-load approach/descend/takeoff overhead on top
                           // of cruise ferry (pilot spirals down & climbs out)
    pub arc_ft_per_turn: f64, // headland turn arc length OUTSIDE the field
                           // polygon billed per turn. Measured median end-to-end
                           // arc 198 ft minus one swath (50 ft) already counted
                           // dead -> ~150 ft net.
    pub ferry_ft: f64,     // one-way trailer -> field
    pub ferry_fps: f64,    // measured 54 mph (incl. climb-out/descent/headland)
    pub t_setup_s: f64,    // flat per-job overhead
    pub t_survey_s: f64,   // field survey, per polygon cluster
    pub survey_link_ft: f64, // polygons closer than this share a survey

    pub dead_fps: f64,       // boom-off transit speed (60 mph)
    pub gap_break_turns: f64, // extra turns a gap turn-around costs per pass
}

impl Default for AircraftParams {
    fn default() -> Self {
        AircraftParams {
            swath_ft: 50.0,
            speed_fps: 83.3,
            t_turn_s: 10.0,
            tank_gal: 75.0,    // bumped 65->75: median good events
                              // had -1 extra load (planner over-predicted)
            gpa: 2.0,
            t_load_s: 100.0,
            t_cycle_s: 10.0,
            arc_ft_per_turn: 150.0,
            ferry_ft: 2640.0,
            ferry_fps: 80.0,
            t_setup_s: 120.0,
            t_survey_s: 180.0,
            survey_link_ft: 1500.0,
            dead_fps: 88.0,
            gap_break_turns: 2.0,
        }
    }
}

impl AircraftParams {
    /// Fly through a dead gap if shorter than this; else turn around.
    /// Breakeven: gap/dead_fps (dead air) vs gap_break_turns * t_turn.
    pub fn gap_threshold_ft(&self) -> f64 {
        self.dead_fps * self.t_turn_s * self.gap_break_turns
    }
}
