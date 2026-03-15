// DedupUMI - Deduplicate FASTQ sequences using UMIs
//
// Rust port of deduplicate_UMI.pl version 1.4
// Original: a.bossers@uu.nl / alex.bossers@wur.nl
// https://github.com/alxelerator/dedupUMI
//
// It no longer uses the system  zlib or zcat/gzip functions but the rust flate2
// I also enabled the zlib-ng SIMD-optimized backend which can be up to 50% faster then zlib by itself.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process;

use chrono::Local;
use clap::Parser;
use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;

const VERSION: &str = "1.4r";
const VERSION_DATE: &str = "2026-03-07";

#[derive(Parser, Debug)]
#[command(name = "deduplicate_UMI", disable_help_flag = true, disable_version_flag = true)]
struct Args {
    #[arg(long = "input-fastq1")]
    input_fastq1: Option<String>,

    #[arg(long = "input-fastq2")]
    input_fastq2: Option<String>,

    #[arg(long = "input-fastq3")]
    input_fastq3: Option<String>,

    #[arg(long = "output-fastq1")]
    output_fastq1: Option<String>,

    #[arg(long = "output-fastq2")]
    output_fastq2: Option<String>,

    #[arg(long = "output-fastq3")]
    output_fastq3: Option<String>,

    #[arg(long = "output-counts")]
    output_counts: Option<String>,

    #[arg(long = "noUMI")]
    no_umi: bool,

    #[arg(long = "verbose")]
    verbose: bool,

    #[arg(long = "help")]
    help: bool,

    #[arg(long = "version")]
    version: bool,
}

fn open_input(path: &str) -> Box<dyn BufRead> {
    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Input file error: {}\n{}", path, e);
        process::exit(1);
    });
    if path.ends_with(".gz") || path.ends_with(".gzip") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

fn open_output(path: &str) -> Box<dyn Write> {
    let file = File::create(path).unwrap_or_else(|e| {
        eprintln!("Output file error: {}\n{}", path, e);
        process::exit(1);
    });
    if path.ends_with(".gz") || path.ends_with(".gzip") {
        Box::new(GzEncoder::new(file, Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    }
}

fn read_line(reader: &mut Box<dyn BufRead>) -> Option<String> {
    let mut buf = String::new();
    match reader.read_line(&mut buf) {
        Ok(0) => None,
        Ok(_) => {
            if buf.ends_with('\n') {
                buf.pop();
                if buf.ends_with('\r') {
                    buf.pop();
                }
            }
            Some(buf)
        }
        Err(e) => {
            eprintln!("Read error: {}", e);
            process::exit(1);
        }
    }
}

fn sum_quality(s: &str) -> u64 {
    s.bytes().map(|b| b as u64).sum()
}

/// Extract UMI from Illumina FASTQ header.
/// Matches Perl regex: ^@.+:[0-9]+:[0-9]+:([ATGCNatgcn]+) [0-9]+:.+
fn extract_umi(header: &str) -> Option<String> {
    if !header.starts_with('@') {
        return None;
    }

    let space_pos = header.find(' ')?;
    let before_space = &header[..space_pos];
    let after_space = &header[space_pos + 1..];

    // Split by ':' to find the last three fields: x_coord, y_coord, UMI
    let parts: Vec<&str> = before_space.split(':').collect();
    if parts.len() < 4 {
        return None;
    }

    let umi = parts[parts.len() - 1];
    let y_coord = parts[parts.len() - 2];
    let x_coord = parts[parts.len() - 3];

    // x and y coordinates must be digit-only strings
    if y_coord.is_empty() || !y_coord.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    if x_coord.is_empty() || !x_coord.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }

    // UMI: non-empty, only ATGCN characters
    if umi.is_empty()
        || !umi
            .chars()
            .all(|c| matches!(c, 'A' | 'T' | 'G' | 'C' | 'N' | 'a' | 't' | 'g' | 'c' | 'n'))
    {
        return None;
    }

    // after_space must match [0-9]+:.+
    let colon_pos = after_space.find(':')?;
    let read_num = &after_space[..colon_pos];
    if read_num.is_empty() || !read_num.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    if after_space[colon_pos + 1..].is_empty() {
        return None;
    }

    Some(umi.to_string())
}

fn print_help() {
    println!("DedupUMI version {} ({})", VERSION, VERSION_DATE);
    println!("https://github.com/alxelerator/dedupUMI");
    println!("a.bossers@uu.nl / alex.bossers@wur.nl");
    println!();
    println!("Usage:");
    println!("  deduplicate_UMI \\");
    println!("      --input-fastq1 <R1.fastq[.gz]> \\");
    println!("      --input-fastq2 <R2.fastq[.gz]> \\");
    println!("      --output-fastq1 <R1_out.fastq[.gz]> \\");
    println!("      --output-fastq2 <R2_out.fastq[.gz]> \\");
    println!("      [--input-fastq3 <UMI.fastq[.gz]>] \\");
    println!("      [--output-fastq3 <UMI_out.fastq[.gz]>] \\");
    println!("      [--output_counts <counts.tab>] \\");
    println!("      [--noUMI] \\");
    println!("      [--verbose] ");
    println!();
    println!("Required arguments:");
    println!("  --input-fastq1 <file>     FASTQ read1 input file (plain or gz)");
    println!("  --input-fastq2 <file>     FASTQ read2 input file (plain or gz)");
    println!("  --output-fastq1 <file>    FASTQ read1 output file (plain or gz)");
    println!("  --output-fastq2 <file>    FASTQ read2 output file (plain or gz)");
    println!();
    println!("Optional arguments:");
    println!("  --input-fastq3 <file>     UMI FASTQ file (R3 system). If omitted, UMI is");
    println!("                            extracted from the read header");
    println!("  --output-fastq3 <file>    Output FASTQ for R3 system UMIs (required if input-fastq3 given)");
    println!("  --output_counts <file>    Append input/output read counts to output tabular file");
    println!("  --noUMI                   Two-file mode without UMI. Deduplicate on R1+R2 only (see README)");
    println!("  --verbose                 Print sequences to STDOUT (debug output)");
    println!("  --help                    Show this help message");
    println!();
    println!("Notes:");
    println!("  - UMI location (header or separate FASTQ file) is detected automatically unless --noUMI is used");
    println!("  - If using --noUMI we just filter exact duplicates. This most likely filters to harsh since");
    println!("    exact duplicates can occurr in natural good quality datasets. Especially at high sequencing depths.");
    println!("  - Gzip I/O is now handled natively via rust::flate2 (zlib-ng backend).");
    println!("    No external gzip/zcat/pigz dependency required anymore from version 1.4r and up.");
    println!();
    println!("Examples:");
    println!("  UMI in header:");
    println!("    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz \\");
    println!("                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz");
    println!("  R3 UMI system:");
    println!("    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --input-fastq3 s1_R3.fastq.gz \\");
    println!("                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz --output-fastq3 s1_R3.dedup.fastq.gz");
    println!("  No UMI:");
    println!("    deduplicate_UMI --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --noUMI \\");
    println!("                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz");
    println!();
}

fn main() {
    let Args {
        input_fastq1,
        input_fastq2,
        input_fastq3,
        output_fastq1,
        output_fastq2,
        output_fastq3,
        output_counts,
        no_umi,
        verbose,
        help,
        version,
    } = Args::parse();

    if help {
        print_help();
        process::exit(0);
    }

    if version {
        println!(
            "DedupUMI version {} ({})\nhttps://github.com/alxelerator/dedupUMI\nBy a.bossers@uu.nl / alex.bossers@wur.nl",
            VERSION, VERSION_DATE
        );
        process::exit(0);
    }

    // Validate required arguments
    if input_fastq1.is_none()
        || input_fastq2.is_none()
        || output_fastq1.is_none()
        || output_fastq2.is_none()
    {
        eprintln!("ERROR: One of the required arguments --input-fastq1, --input-fastq2, --input-fastq3 or --output-fastq1, --output-fastq2, --output-fastq3 is missing!");
        println!("DedupUMI version {} ({})", VERSION, VERSION_DATE);
        println!("Use: deduplicate_UMI --help for options");
        process::exit(1);
    }

    // UMI compatibility checks
    if no_umi && input_fastq3.is_some() {
        eprintln!("ERROR: --noUMI cannot be combined with --input-fastq3");
        process::exit(1);
    }
    if no_umi && output_fastq3.is_some() {
        eprintln!("ERROR: --noUMI cannot be combined with --output-fastq3");
        process::exit(1);
    }

    // Open counts output file (append if exists, create with header if new)
    let mut counts_writer: Option<BufWriter<File>> = if let Some(ref path) = output_counts {
        if Path::new(path).exists() {
            println!("dedupUMI: Appending read count data to existing file!");
            let f = OpenOptions::new()
                .append(true)
                .open(path)
                .unwrap_or_else(|e| {
                    eprintln!("Output file error: {}\n{}", path, e);
                    process::exit(1);
                });
            Some(BufWriter::new(f))
        } else {
            let mut f = File::create(path).unwrap_or_else(|e| {
                eprintln!("Output file error: {}\n{}", path, e);
                process::exit(1);
            });
            writeln!(f, "Sample\tReads_in\tReads_out").unwrap_or_else(|e| {
                eprintln!("Output file error: {}\n{}", path, e);
                process::exit(1);
            });
            Some(BufWriter::new(f))
        }
    } else {
        None
    };

    // Determine processing mode
    let umi_header: bool;
    if input_fastq3.is_some() {
        if output_fastq3.is_none() {
            eprintln!("DeDupe ERROR: R3-system (R1 R2 R3) detected, but output file R3 is not defined! Define by using --output-fastq3");
            process::exit(1);
        }
        println!("DeDupe: UMI sequences processing R3-system (R1 R2 R3)");
        umi_header = false;
    } else if no_umi {
        println!("DeDupe: Processing noUMI mode (R1 R2 exact duplicate read-pairs)");
        umi_header = false;
    } else {
        println!("DeDupe: UMI sequences processing UMI in the R1 R2 headers!");
        umi_header = true;
    }

    // Build input/output file lists
    let input_fastq1 = input_fastq1.unwrap();
    let mut infiles: Vec<String> = vec![input_fastq1.clone(), input_fastq2.unwrap()];
    let mut outfiles: Vec<String> = vec![output_fastq1.unwrap(), output_fastq2.unwrap()];

    if let Some(f3) = input_fastq3 {
        infiles.push(f3);
        outfiles.push(output_fastq3.unwrap());
    }

    let nfiles = infiles.len();

    let mut readers: Vec<Box<dyn BufRead>> = infiles.iter().map(|f| open_input(f)).collect();
    let mut writers: Vec<Box<dyn Write>> = outfiles.iter().map(|f| open_output(f)).collect();

    // Main read loop: store best (highest quality) record per unique key
    let mut reads_in: u64 = 0;
    // HashMap key   -> (head_joined, qual_joined)
    // head_joined   = headers[0] + "#alx#" + headers[1] [+ "#alx#" + headers[2]]
    // qual_joined   = quality[0] + " " + quality[1] [+ " " + quality[2]]
    let mut uniqseq: HashMap<String, (String, String)> = HashMap::new();

    loop {
        // Read header line from file 0
        let header0 = match read_line(&mut readers[0]) {
            Some(h) => h,
            None => break,
        };
        reads_in += 1;

        // Read header lines from remaining files
        let mut headers = vec![header0];
        for i in 1..nfiles {
            match read_line(&mut readers[i]) {
                Some(h) => headers.push(h),
                None => {
                    eprintln!(
                        "ERROR: Paired FASTQ ended early / files out of sync (file {} missing header while R1 still has data)",
                        i + 1
                    );
                    process::exit(1);
                }
            }
        }

        // Empty header signals trailing newline at EOF
        if headers[0].is_empty() {
            eprintln!("WARNING: Empty header line in R1 (extra newline at EOF?)");
            break;
        }

        // Read sequence lines
        let mut sequences: Vec<String> = Vec::with_capacity(nfiles);
        for i in 0..nfiles {
            match read_line(&mut readers[i]) {
                Some(s) => sequences.push(s),
                None => {
                    eprintln!(
                        "ERROR: FASTQ ended early / missing sequence line in file {}",
                        i + 1
                    );
                    process::exit(1);
                }
            }
        }

        // Extract UMI from header (UMI-in-header mode only)
        let umi = if umi_header {
            match extract_umi(&headers[0]) {
                Some(u) => u,
                None => {
                    eprintln!(
                        "  ERROR: UMI not found in fastq header or fastq header format has changed!?\n    Header: '{}'",
                        headers[0]
                    );
                    eprintln!("  Consider using the option --noUMI if your sequence library did not contain UMI sequences!");
                    process::exit(1);
                }
            }
        } else {
            String::new()
        };

        // Skip '+' lines
        for i in 0..nfiles {
            match read_line(&mut readers[i]) {
                Some(_) => {}
                None => {
                    eprintln!(
                        "ERROR: FASTQ ended early / missing '+' line in file {}",
                        i + 1
                    );
                    process::exit(1);
                }
            }
        }

        // Read quality lines
        let mut qualities: Vec<String> = Vec::with_capacity(nfiles);
        for i in 0..nfiles {
            match read_line(&mut readers[i]) {
                Some(q) => qualities.push(q),
                None => {
                    eprintln!(
                        "ERROR: FASTQ ended early / missing quality line in file {}",
                        i + 1
                    );
                    process::exit(1);
                }
            }
        }

        // Build deduplication key: sequences joined by space (+ UMI if applicable)
        let key = if umi_header {
            format!("{} {} {}", sequences[0], sequences[1], umi)
        } else if nfiles == 3 {
            format!("{} {} {}", sequences[0], sequences[1], sequences[2])
        } else {
            format!("{} {}", sequences[0], sequences[1])
        };

        // If duplicate exists, keep the record with higher total quality score
        if let Some((_, stored_qual)) = uniqseq.get(&key) {
            let org_quals: Vec<&str> = stored_qual.split(' ').collect();

            let (total, org_total) = if umi_header || nfiles != 3 {
                // UMI-in-header or noUMI: compare R1 + R2 quality
                let t = sum_quality(&qualities[0]) + sum_quality(&qualities[1]);
                let ot = sum_quality(org_quals[0]) + sum_quality(org_quals[1]);
                (t, ot)
            } else {
                // R3 system: compare R1 + R3 quality (UMI is in R3 file)
                let t = sum_quality(&qualities[0]) + sum_quality(&qualities[2]);
                let ot = sum_quality(org_quals[0]) + sum_quality(org_quals[2]);
                (t, ot)
            };

            if total <= org_total {
                continue; // keep existing record
            }
        }

        // Store (or replace) the record
        uniqseq.insert(key, (headers.join("#alx#"), qualities.join(" ")));
    }

    println!("DeDupe: Read    : {} sequences", reads_in);

    // Write all unique records to output files
    let mut seqcount_w: u64 = 0;

    for (seqkey, (head_joined, qual_joined)) in &uniqseq {
        if verbose {
            println!("  {}", seqkey);
        }

        let seqs: Vec<&str> = seqkey.split(' ').collect();
        let hdrs: Vec<&str> = head_joined.split("#alx#").collect();
        let quals: Vec<&str> = qual_joined.split(' ').collect();

        for i in 0..nfiles {
            let w = &mut writers[i];
            writeln!(w, "{}", hdrs[i]).unwrap_or_else(|e| {
                eprintln!("Write error: {}", e);
                process::exit(1);
            });
            writeln!(w, "{}", seqs[i]).unwrap_or_else(|e| {
                eprintln!("Write error: {}", e);
                process::exit(1);
            });
            writeln!(w, "+").unwrap_or_else(|e| {
                eprintln!("Write error: {}", e);
                process::exit(1);
            });
            writeln!(w, "{}", quals[i]).unwrap_or_else(|e| {
                eprintln!("Write error: {}", e);
                process::exit(1);
            });
        }

        seqcount_w += 1;
    }

    println!("Dedupe: Written : {} sequences", seqcount_w);
    println!(
        "DeDupe: Finished {}\n",
        Local::now().format("%m/%d/%Y %H:%M:%S")
    );

    // Append counts to tabular file if requested
    if let Some(ref mut cw) = counts_writer {
        writeln!(cw, "{}\t{}\t{}", input_fastq1, reads_in, seqcount_w).unwrap_or_else(|e| {
            eprintln!("Counts write error: {}", e);
            process::exit(1);
        });
    }
}
