use bio_rascal::gff::parse_gff_line;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(c: &mut Criterion) {
    let test_line = include_str!("../src/test-data/test-line.gff")
        .trim()
        .to_string();
    c.bench_function("test skip", |b| {
        b.iter(|| parse_gff_line(black_box(&test_line), false))
    });
    c.bench_function("test include", |b| {
        b.iter(|| parse_gff_line(black_box(&test_line), true))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
