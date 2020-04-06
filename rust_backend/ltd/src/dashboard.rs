use integrand::IntegrandStatistics;
use observables::EventInfo;
use std::io;
use std::sync::mpsc::{channel, Sender};
use std::thread;
use std::time::Duration;
use tui::layout::{Constraint, Corner, Direction, Layout};
use tui::style::{Color, Modifier, Style};
use tui::widgets::{Axis, Block, Borders, Chart, Dataset, List, Marker, Text};
use tui::{backend::TermionBackend, Terminal};

pub type StatusUpdateSender = Sender<StatusUpdate>;

#[derive(Clone, Debug)]
pub enum StatusUpdate {
    NewPoint(usize, f64, f64, f64, f64, f64, f64),
    IntegratorUpdate(String),
    Statistics(IntegrandStatistics),
    EventInfo(EventInfo),
    Message(String),
}

pub struct Dashboard {
    pub status_update_sender: Sender<StatusUpdate>,
}

impl Dashboard {
    pub fn new(full: bool) -> Dashboard {
        if full {
            Dashboard::full_dashboard()
        } else {
            Dashboard::minimal_dashboard()
        }
    }

    pub fn minimal_dashboard() -> Dashboard {
        let (log_sender, log_receiver) = channel();

        thread::spawn(move || loop {
            while let Ok(x) = log_receiver.recv() {
                match x {
                    StatusUpdate::Message(m) => println!("{}", m),
                    StatusUpdate::IntegratorUpdate(_) => {}
                    StatusUpdate::NewPoint(..) => {}
                    StatusUpdate::Statistics(s) => {
                        println!("{:?}", s);
                    }
                    StatusUpdate::EventInfo(s) => {
                        if s.accepted_event_counter != 0 {
                            println!("{:?}", s);
                        }
                    }
                }
            }
        });

        Dashboard {
            status_update_sender: log_sender,
        }
    }

    fn full_dashboard() -> Dashboard {
        let (log_sender, log_receiver) = channel();

        let stdout = io::stdout(); //.into_raw_mode()?;
        let backend = TermionBackend::new(stdout);
        let mut terminal = Terminal::new(backend).unwrap();
        terminal.hide_cursor().unwrap(); // this does not get reset after the program closes

        let mut re_data = vec![];
        let mut im_data = vec![];

        let mut log_messages = vec![];
        let mut integrator_update_messages = vec![];

        let mut integrand_statistics = IntegrandStatistics::new();
        let mut event_info = EventInfo::default();

        terminal.clear().unwrap();
        thread::spawn(move || loop {
            while let Ok(x) = log_receiver.try_recv() {
                match x {
                    StatusUpdate::Message(m) => log_messages.push(Text::raw(m)),
                    StatusUpdate::IntegratorUpdate(m) => {
                        integrator_update_messages.push(Text::raw(m))
                    }
                    StatusUpdate::NewPoint(iter, re, re_err, re_chi, im, im_err, im_chi) => {
                        re_data.push((iter as f64, re));
                        im_data.push((iter as f64, im));

                        integrator_update_messages.push(Text::raw(format!("Iteration {}", iter)));
                        integrator_update_messages.push(Text::raw(format!(
                            " re: {:+.8e} +- {:+.8e} chisq {:.2}",
                            re, re_err, re_chi
                        )));
                        integrator_update_messages.push(Text::raw(format!(
                            " im: {:+.8e} +- {:+.8e} chisq {:.2}",
                            im, im_err, im_chi
                        )));
                    }
                    StatusUpdate::EventInfo(e) => {
                        event_info = e;
                    }
                    StatusUpdate::Statistics(s) => {
                        integrand_statistics = s;
                    }
                }
            }
            terminal
                .draw(|mut f| {
                    let vert_chunks = Layout::default()
                        .direction(Direction::Vertical)
                        .constraints(
                            [Constraint::Percentage(80), Constraint::Percentage(20)].as_ref(),
                        )
                        .split(f.size());
                    let chunks_hor = Layout::default()
                        .direction(Direction::Horizontal)
                        .constraints(
                            [
                                Constraint::Percentage(25),
                                Constraint::Percentage(25),
                                Constraint::Percentage(50),
                            ]
                            .as_ref(),
                        )
                        .split(vert_chunks[0]);

                    let mut events_list =
                        List::new(integrator_update_messages.iter().rev().cloned())
                            .block(
                                Block::default()
                                    .borders(Borders::ALL)
                                    .title("Integrator updates"),
                            )
                            .start_corner(Corner::BottomLeft);
                    f.render(&mut events_list, chunks_hor[0]);

                    let stats = vec![
                        Text::raw(format!(
                            "Total points: {}",
                            integrand_statistics.total_samples
                        )),
                        Text::styled(
                            format!(
                                "Unstable points: {} ({:.2}%)",
                                integrand_statistics.unstable_point_count,
                                integrand_statistics.unstable_point_count as f64
                                    / integrand_statistics.total_samples as f64
                                    * 100.
                            ),
                            Style::default().fg(Color::Yellow),
                        ),
                        Text::styled(
                            format!(
                                "Unstable quad points: {} ({:.2}%)",
                                integrand_statistics.unstable_f128_point_count,
                                integrand_statistics.unstable_f128_point_count as f64
                                    / integrand_statistics.total_samples as f64
                                    * 100.
                            ),
                            Style::default().fg(Color::Red),
                        ),
                        Text::raw(format!(
                            "Evaluation time per sample: {:.2}µs",
                            integrand_statistics.total_sample_time
                                / integrand_statistics.total_samples as f64
                        )),
                        Text::raw(format!("Running max: {:e}", integrand_statistics.running_max)),
                        Text::raw(format!(
                            "Accepted events: {}",
                            event_info.accepted_event_counter
                        )),
                        Text::raw(format!(
                            "Rejected events: {} ({:.2}%)",
                            event_info.rejected_event_counter,
                            event_info.rejected_event_counter as f64
                                / (event_info.rejected_event_counter as f64
                                    + event_info.accepted_event_counter as f64)
                                * 100.
                        )),
                    ];

                    let mut stats_list = List::new(stats.into_iter())
                        .block(
                            Block::default()
                                .borders(Borders::ALL)
                                .title("Integrator statistics"),
                        )
                        .start_corner(Corner::TopLeft);
                    f.render(&mut stats_list, chunks_hor[1]);

                    let datasets = [
                        Dataset::default()
                            .name("re")
                            .marker(Marker::Dot)
                            .style(Style::default().fg(Color::Cyan))
                            .data(re_data.as_slice()),
                        Dataset::default()
                            .name("im")
                            .marker(Marker::Braille)
                            .style(Style::default().fg(Color::Yellow))
                            .data(im_data.as_slice()),
                    ];
                    let x_labels = [
                        "0".to_owned(),
                        (re_data.len() / 2).to_string(),
                        re_data.len().to_string(),
                    ];

                    let bounds = [
                        re_data
                            .iter()
                            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                            .map(|x| x.1)
                            .unwrap_or(0.),
                        re_data
                            .iter()
                            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                            .map(|x| x.1)
                            .unwrap_or(0.),
                    ];
                    let y_labels = [
                        format!("{:.2e}", bounds[0]),
                        format!("{:.2e}", (bounds[0] + bounds[1]) / 2.),
                        format!("{:.2e}", bounds[1]),
                    ];

                    let mut chart = Chart::default()
                        .block(
                            Block::default()
                                .title("Integrand")
                                .title_style(
                                    Style::default().fg(Color::Cyan).modifier(Modifier::BOLD),
                                )
                                .borders(Borders::NONE),
                        )
                        .x_axis(
                            Axis::default()
                                .title("Iteration")
                                .style(Style::default().fg(Color::Gray))
                                .labels_style(Style::default().modifier(Modifier::ITALIC))
                                .bounds([0., re_data.len() as f64])
                                .labels(&x_labels),
                        )
                        .y_axis(
                            Axis::default()
                                .labels_style(Style::default().modifier(Modifier::ITALIC))
                                .bounds(bounds)
                                .labels(&y_labels),
                        )
                        .datasets(&datasets);
                    f.render(&mut chart, chunks_hor[2]);

                    let block = Block::default()
                        .borders(Borders::ALL)
                        .title("Log")
                        .title_style(Style::default().fg(Color::Magenta).modifier(Modifier::BOLD));
                    let mut log = List::new(log_messages.iter().rev().cloned())
                        .block(block)
                        .start_corner(Corner::BottomLeft);
                    f.render(&mut log, vert_chunks[1]);
                })
                .unwrap();
            thread::sleep(Duration::from_millis(200));
        });

        Dashboard {
            status_update_sender: log_sender,
        }
    }
}
