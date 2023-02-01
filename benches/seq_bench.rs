use bio_rascal::sequence::{complement, reverse_comp, universal_code, UNIVERSAL_CODE};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(c: &mut Criterion) {
    let original = b"ACTGATATATGCGCGCATCTC";
    c.bench_function("rev/comp", |b| {
        b.iter(|| complement(black_box(original.iter().rev())))
    });

    c.bench_function("reverse_comp", |b| {
        b.iter(|| reverse_comp(black_box(original)))
    });

    let codons: Vec<Vec<u8>> = UNIVERSAL_CODE
        .keys()
        .map(|k| k.as_bytes().to_vec())
        .collect();
    c.bench_function("test &[u8] to function", |b| {
        b.iter(|| universal_code(black_box(&codons[0])))
    });
    c.bench_function("test &[u8] to phf map", |b| {
        b.iter(|| UNIVERSAL_CODE.get(black_box(std::str::from_utf8(&codons[0]).unwrap())))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
