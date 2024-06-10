use std::fmt::Display;

use mathru::algebra::linear::{
    matrix::{General, Solve},
    vector::Vector,
};

use crate::{Point, TNumber};

pub trait Function: Display {
    fn compute(&self, x: TNumber) -> TNumber;
}

pub trait MinimizedFunction {
    fn new_minimized(points: &Vec<Point>) -> Self;
}

pub fn create_approximations(points: &Vec<Point>) -> Vec<Box<dyn Function>> {
    vec![
        Box::new(Linear::new_minimized(points)),
        Box::new(Quadratic::new_minimized(points)),
        Box::new(Cubic::new_minimized(points)),
        Box::new(Exponent::new_minimized(points)),
        Box::new(Logrithm::new_minimized(points)),
        Box::new(Power::new_minimized(points)),
    ]
}

pub struct Linear {
    /// Multiplier
    a: TNumber,
    /// Addition
    b: TNumber,
}

// special thanks to Lannee for implementation of all minimization rutines

impl Display for Linear {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Linear")?;
        writeln!(f, "{}*x + {}", self.a, self.b)
    }
}

impl Function for Linear {
    fn compute(&self, x: TNumber) -> TNumber {
        self.a * x + self.b
    }
}

impl MinimizedFunction for Linear {
    fn new_minimized(points: &Vec<Point>) -> Linear {
        let (sx, sxx, sy, sxy) = points
            .iter()
            .fold((0., 0., 0., 0.), |(sx, sxx, sy, sxy), Point { x, y }| {
                (sx + x, sxx + x.powi(2), sy + y, sxy + x * y)
            });

        let n = points.len() as f64;
        let a = (sxy * n - sx * sy) / (sxx * n - sx.powi(2));
        let b = (sxx * sy - sx * sxy) / (sxx * n - sx.powi(2));

        Linear { a, b }
    }
}

pub struct Quadratic {
    a0: TNumber,
    a1: TNumber,
    a2: TNumber,
}

impl Display for Quadratic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Quadratic")?;
        writeln!(f, "{} + {}*x + {}*x^2", self.a0, self.a1, self.a2)
    }
}

impl Function for Quadratic {
    fn compute(&self, x: TNumber) -> TNumber {
        self.a0 + self.a1 * x + self.a2 * x.powi(2)
    }
}

impl MinimizedFunction for Quadratic {
    fn new_minimized(points: &Vec<Point>) -> Self {
        let mut matrix = General::<f64>::zero(3, 3);
        let mut vector = Vector::<f64>::zero(3);

        points.iter().for_each(|&Point { x, y }| {
            matrix[[0, 0]] += 1.;
            matrix[[0, 1]] += x;
            matrix[[0, 2]] += x.powi(2);
            matrix[[1, 0]] += x;
            matrix[[1, 1]] += x.powi(2);
            matrix[[1, 2]] += x.powi(3);
            matrix[[2, 0]] += x.powi(2);
            matrix[[2, 1]] += x.powi(3);
            matrix[[2, 2]] += x.powi(4);

            vector[0] += y;
            vector[1] += x * y;
            vector[2] += x * x * y;
        });

        let coeffs = matrix.solve(&vector).unwrap();
        let a0 = coeffs[0];
        let a1 = coeffs[1];
        let a2 = coeffs[2];

        Quadratic { a0, a1, a2 }
    }
}

pub struct Cubic {
    a0: TNumber,
    a1: TNumber,
    a2: TNumber,
    a3: TNumber,
}

impl Display for Cubic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Cubic")?;
        writeln!(
            f,
            "{} + {}*x + {}*x^2 + {}*x^3",
            self.a0, self.a1, self.a2, self.a3
        )
    }
}

impl Function for Cubic {
    fn compute(&self, x: TNumber) -> TNumber {
        self.a0 + self.a1 * x + self.a2 * x.powi(2) + self.a3 * x.powi(3)
    }
}

impl MinimizedFunction for Cubic {
    fn new_minimized(points: &Vec<Point>) -> Self {
        let mut matrix = General::<f64>::zero(4, 4);
        let mut vector = Vector::<f64>::zero(4);

        points.iter().for_each(|&Point { x, y }| {
            matrix[[0, 0]] += 1.;
            matrix[[0, 1]] += x;
            matrix[[0, 2]] += x.powi(2);
            matrix[[0, 3]] += x.powi(3);
            matrix[[1, 0]] += x;
            matrix[[1, 1]] += x.powi(2);
            matrix[[1, 2]] += x.powi(3);
            matrix[[1, 3]] += x.powi(4);
            matrix[[2, 0]] += x.powi(2);
            matrix[[2, 1]] += x.powi(3);
            matrix[[2, 2]] += x.powi(4);
            matrix[[2, 3]] += x.powi(5);
            matrix[[3, 0]] += x.powi(3);
            matrix[[3, 1]] += x.powi(4);
            matrix[[3, 2]] += x.powi(5);
            matrix[[3, 3]] += x.powi(6);

            vector[0] += y;
            vector[1] += x * y;
            vector[2] += x.powi(2) * y;
            vector[3] += x.powi(3) * y;
        });

        let coeffs = matrix.solve(&vector).unwrap();
        let a0 = coeffs[0];
        let a1 = coeffs[1];
        let a2 = coeffs[2];
        let a3 = coeffs[3];

        Cubic { a0, a1, a2, a3 }
    }
}

pub struct Exponent {
    a0: TNumber,
    a1: TNumber,
}

impl Display for Exponent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Exponent")?;
        writeln!(f, "e^({}*x + {})", self.a0, self.a1)
    }
}

impl Function for Exponent {
    fn compute(&self, x: TNumber) -> TNumber {
        (self.a0 * x + self.a1).exp()
    }
}

impl MinimizedFunction for Exponent {
    fn new_minimized(points: &Vec<Point>) -> Self {
        let points: Vec<_> = points
            .iter()
            .map(|Point { x, y }| Point { x: *x, y: y.ln() })
            .collect();

        let Linear { a: a0, b: a1 } = Linear::new_minimized(&points);

        Exponent { a0, a1 }
    }
}

pub struct Logrithm {
    a0: TNumber,
    a1: TNumber,
}

impl Display for Logrithm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Natural Logarithm")?;
        writeln!(f, "{} * ln(x) + {}", self.a0, self.a1)
    }
}

impl Function for Logrithm {
    fn compute(&self, x: TNumber) -> TNumber {
        self.a0 * x.ln() + self.a1
    }
}

impl MinimizedFunction for Logrithm {
    fn new_minimized(points: &Vec<Point>) -> Self {
        let points_mapped: Vec<_> = points
            .iter()
            .map(|&Point { x, y }| Point { x: x.ln(), y })
            .collect();

        let Linear { a: a0, b: a1 } = Linear::new_minimized(&points_mapped);

        Logrithm { a0, a1 }
    }
}

pub struct Power {
    a0: TNumber,
    a1: TNumber,
}

impl Display for Power {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Power")?;
        writeln!(f, "{}*e*x^{}", self.a0, self.a1)
    }
}

impl Function for Power {
    fn compute(&self, x: TNumber) -> TNumber {
        self.a0.exp() * x.powf(self.a1)
    }
}

impl MinimizedFunction for Power {
    fn new_minimized(points: &Vec<Point>) -> Self {
        let points_mapped: Vec<_> = points
            .iter()
            .map(|Point { x, y }| Point {
                x: x.ln(),
                y: y.ln(),
            })
            .collect();

        let Linear { a: a0, b: a1 } = Linear::new_minimized(&points_mapped);

        Power { a0, a1 }
    }
}
