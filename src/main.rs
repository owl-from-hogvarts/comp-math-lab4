use std::{
    error::Error,
    ops::{Deref, Range},
};

use cli_table::format::{Border, HorizontalLine, Separator, VerticalLine};
use methods::Function;
use serde::Deserialize;

use crate::methods::create_approximations;

mod methods;

type TNumber = f64;

#[derive(Clone, Copy, Debug, Deserialize)]
struct Point {
    pub x: TNumber,
    pub y: TNumber,
}

fn main() {
    match start() {
        Ok(_) => (),
        Err(error) => eprintln!("{}", error),
    }
}

// ln(0) = -inf
// so we need to account for point (0, y)
const APPROX_ZERO: TNumber = 0.000001;

fn start() -> Result<(), Box<dyn Error>> {
    let points: Vec<Point> = input_points()?
        .iter()
        .map(|&Point { x, y }| Point {
            x: if x == 0. { APPROX_ZERO } else { x },
            y: if y == 0. { APPROX_ZERO } else { y },
        })
        .collect();

    // compute minimal for each funciton
    let all_approximations = create_approximations(&points);

    // compute total deviation for each
    let approximated_points: Vec<_> = all_approximations
        .iter()
        .map(|function| compute_deviation(&points, function.deref()))
        .collect();
    let standard_deviations: Vec<f64> = approximated_points
        .iter()
        .map(|deviations| {
            deviations
                .iter()
                .map(|(_, _, epsilon)| epsilon.powi(2))
                .sum::<f64>()
        })
        .map(|epsilon_sum| (epsilon_sum / points.len() as f64).sqrt())
        .collect();

    let best_approximation = standard_deviations
        .iter()
        .zip(all_approximations)
        .enumerate()
        .min_by(move |(_, a), (_, b)| a.0.total_cmp(b.0))
        .expect("At least one approximation present");

    println!("{}", best_approximation.1 .1);
    println!("Standard deviation is: {:.5}", best_approximation.1 .0);
    print_points(approximated_points.get(best_approximation.0).expect(
        "amount of approximation arrays should match with number of approximation functions",
    ))?;

    plot(&points, best_approximation.1 .1.deref())
}

fn input_points() -> Result<Vec<Point>, serde_json::Error> {
    serde_json::from_reader(std::io::stdin())
}

fn compute_deviation(
    points: &Vec<Point>,
    function: &dyn Function,
) -> Vec<(Point, TNumber, TNumber)> {
    points
        .iter()
        .map(|&point| {
            let phi = function.compute(point.x);
            let epsilon = phi - point.y;
            (point, phi, epsilon)
        })
        .collect()
}

fn print_points(points: &Vec<(Point, f64, f64)>) -> Result<(), Box<dyn std::error::Error>> {
    use cli_table::Table;
    let table = points
        .iter()
        .enumerate()
        .map(|(index, point)| {
            vec![
                (index + 1).to_string(),
                format!("{:.4}", point.0.x),
                format!("{:.4}", point.0.y),
                format!("{:.4}", point.1),
                format!("{:.4}", point.2),
            ]
        })
        .table()
        .border(
            Border::builder()
                .top(HorizontalLine::new('╭', '╮', '┬', '─'))
                .left(VerticalLine::new('│'))
                .right(VerticalLine::new('│'))
                .bottom(HorizontalLine::new('╰', '╯', '┴', '─'))
                .build(),
        )
        .separator(
            Separator::builder()
                .row(Some(HorizontalLine::new('├', '┤', '┼', '─')))
                .build(),
        )
        .title(["Point number", "X", "Y", "φ(x)", "ε"])
        .display()?;

    println!("{table}");
    return Ok(());
}

fn with_coord_margin(range: Range<f64>, margin_persents: f64) -> Range<f64> {
    let length = range.end - range.start;
    let margin = length * margin_persents;
    (range.start - margin)..(range.end + margin)
}

fn plot(points: &Vec<Point>, function: &dyn Function) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;
    const MARGINS: i32 = 10;
    const COORD_MARGIN_PERSENT: TNumber = 0.05;
    const IMAGE_PATH: &'static str = "./plot.png";
    const POINT_SIZE: i32 = 10;

    println!("Generating image. This may take several seconds");

    let x_range = {
        let min = points
            .iter()
            .min_by(|a, b| a.x.total_cmp(&b.x))
            .expect("At least one point present");
        let max = points
            .iter()
            .max_by(|a, b| a.x.total_cmp(&b.x))
            .expect("At least one point present");

        min.x..max.x
    };

    let y_range = {
        let min = points
            .iter()
            .min_by(|a, b| a.y.total_cmp(&b.y))
            .expect("At least one point present");
        let max = points
            .iter()
            .max_by(|a, b| a.y.total_cmp(&b.y))
            .expect("At least one point present");

        min.y..max.y
    };

    let root = BitMapBackend::new(IMAGE_PATH, (1920, 1080)).into_drawing_area();
    root.fill(&WHITE)?;
    let root = root.margin(MARGINS, MARGINS, MARGINS, MARGINS);

    let mut chart = ChartBuilder::on(&root)
        .margin(MARGINS * 2)
        .x_label_area_size(20)
        .y_label_area_size(40)
        .build_cartesian_2d(
            with_coord_margin(x_range.clone(), COORD_MARGIN_PERSENT),
            with_coord_margin(y_range, COORD_MARGIN_PERSENT),
        )?;

    chart
        .configure_mesh()
        .label_style(("noto sans", 16))
        .x_labels(5)
        .y_labels(5)
        .x_desc("X")
        .y_desc("Y")
        .draw()?;
    chart.draw_series(PointSeries::<_, _, Circle<_, _>, _>::new(
        points.iter().map(|point| (point.x, point.y)),
        POINT_SIZE,
        BLACK.filled(),
    ))?;

    chart.draw_series(LineSeries::new(
        x_range
            .clone()
            .step(0.05)
            .values()
            .chain([x_range.end])
            .map(|x| (x, function.compute(x))),
        GREEN.stroke_width(3),
    ))?;

    root.present()?;

    println!("Image saved at path: {}", IMAGE_PATH);
    Ok(())
}
