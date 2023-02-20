# A bioinformatics library in Rust

**Highly experimental**

This library serves 2 purposes:

1) helps me learn to code in Rust by porting existing Python code
2) use it as a base for new programs I need/want to make

## Why Rust

The vast majority of software that I code for bioinformatics is CLI, runs on a HPC o remote computer. In general I used Python for the base of my research and maintained a [library](https://github.com/frubino/mgkit) for a long time. Python is amazing for a lot of project, especially when it comes to data analysis and the advantages are a lot. However, some things have been a problem for a while on a specific HPC:

- size of distribution, especially using Conda
- startup time, because of the high volume of files and library loaded

These can be mitigated with Rust, especially using a static linking. Moreover, memory and speed requirements for large datasets have required a lot of workarounds. So it makes sense for me to learn a new compiled language and start porting some utilities from the library above.

## Why (or why not) use this library

While I am using this library for some software I am making, it may be useful for other people. Having said that, I am changing a lot inside it as I use in other projects and **it's not stable at the moment**. I think the objective for version 0.2.0 will be to try and stabilise the API along with the other software. In the meantime, it can be useful for someone else to experiment.

## Why not using library *X* for *Y*

One reason is that I am learning a new language and I find it much easier to adapt code I wrote before in another language. In fact, my understanding of Rust became much better when I started using the basics, instead of fighting a particular API. Also, these libraries were a bit difficult to get into when I started but I will probably move to use portions of them when I feel more confident using Rust.

My primary objective is writing a few programs I need and I need a library to reuse code. This allows me to learn more testing, benchmarking and packaging in Rust. In the end, a learning experience and I tried the same with Go as an alternative, but I prefer Rust so far.

## Testing

There are tests for the library and some integration tests. You can run the test with `cargo test` or `cargo nextest run` if [nextest](https://nexte.st/) is installed. The integration test expect a taxonomy file `taxonomy.json.gz` in `tests/test-data/`, which can be created using [taxa-utils](https://github.com/frubino/taxa-utils).

## Benchmarks

There are only a few that I used to test some changes and to learn them. They require [criterion](https://github.com/bheisler/criterion.rs) to run.