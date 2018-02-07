extern crate rand;
extern crate clap;

use clap::{Arg, App};
use std::fmt;
use rand::Rng;
use std::collections::BTreeSet;

fn exit(msg: &str) {
    eprintln!("{}", msg);
    std::process::exit(1);
}

fn is_usize(s: String) -> Result<(), String> {
    match s.parse::<usize>() {
        Ok(u) => Ok(()),
        Err(_) => Err("Ist kein Integer-Wert".to_string())
    }
}

fn main() {
    let matches = App::new("btwns")
        .version("1.0")
        .about("Bearbeitet das Problem MAXIMUM BETWEENNESS mittels evolutionärem Algorithmus")
        .author("Felix Haller <felix.haller@stud.htwk-leipzig.de>")
        .arg(Arg::with_name("mutation")
            .short("m")
            .value_name("MUTATION")
            .help("Mutationsstrategie")
            .possible_values(&["swap", "shift", "invert"])
            .default_value("swap")
            .takes_value(true)
            .required(false))
        .arg(Arg::with_name("elements")
            .short("e")
            .value_name("ELEMENTS")
            .help("Anzahl der Elemente der Permutation")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("100"))
        .arg(Arg::with_name("triples")
            .short("t")
            .value_name("TRIPLES")
            .help("Anzahl der Triples")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("50"))
        .arg(Arg::with_name("generations")
            .short("g")
            .value_name("GENERATIONS")
            .help("Anzahl max. Generationen")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("1000"))
        .arg(Arg::with_name("unchanged")
            .short("u")
            .value_name("UNCHANGED")
            .help("Anzahl der max. aufeinanderfolgenden Generationen ohne Veränderung")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("100"))
        .arg(Arg::with_name("popsize")
            .short("p")
            .long("pop")
            .value_name("POPULATIONSIZE")
            .help("Populationsgröße festlegen")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("100"))
        .get_matches();

    let mutation = matches.value_of("mutation").unwrap();
    let elements = matches.value_of("elements").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let triples = matches.value_of("triples").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let generations = matches.value_of("generations").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let unchanged = matches.value_of("unchanged").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let popsize = matches.value_of("popsize").unwrap().parse::<usize>().expect("Fehler beim Parsen");

    solve_solvable(elements, triples, generations, unchanged, popsize, mutation);
}

struct Individual {
    id : u32,
    generation: u32,
    // Das Genome
    genome: Vec<usize>,
    // Eine umgekehrte Zuordnung um schnell den Index eines Elements zu erfahren
}

    impl fmt::Display for Individual {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "#{}[{}]: {}", self.id, self.generation, format!("{:?}", self.genome))
        }
    }

    impl Individual {
        fn new_random(elements: usize) -> Individual {
            let mut genome = Vec::new();
            for i in 0..elements {
                genome.push(i);
            }
            rand::thread_rng().shuffle(&mut genome);


            Individual{id : 0, generation : 1, genome: genome}
        }

        fn mutate(&mut self) {

            // vertauschende Mutation
            let a = rand::thread_rng().gen_range(0, self.genome.len());
            let b = rand::thread_rng().gen_range(0, self.genome.len());

            self.genome.swap(a, b);


            // invertierende Mutation

            // verschiebende Mutation (gut geeignet laut Skript)

        }
        fn recombine(&self, other : Individual) {
            // ORDNUNGSREKOMBINATION (Buch S.29)
        }
    }

struct BTWNSProblem {
    rules : BTreeSet<(usize,usize,usize)>,
}

    impl fmt::Display for BTWNSProblem {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{:?}", self.rules)
        }
    }

    impl BTWNSProblem {
        fn new_solvable_from_individual(solution: &Individual, elements: usize, triples: usize) -> BTWNSProblem {
            let mut rules = BTreeSet::<(usize,usize,usize)>::new();
            while rules.len() < triples {
                let i_c = rand::thread_rng().gen_range(2, elements);
                let i_b = rand::thread_rng().gen_range(1, i_c);
                let i_a = rand::thread_rng().gen_range(0, i_b);
                let c = solution.genome[i_c];
                let b = solution.genome[i_b];
                let a = solution.genome[i_a];
                rules.insert((a,b,c));
            }
            BTWNSProblem{rules : rules}
        }

        fn rate(&self, ind: &Individual) -> f64 {

            // genome und r_genome erzeugen um später in konstanter Zeit die Position eines
            // Elements zu erfahren
            let mut r_genome = vec![0; ind.genome.len()];
            for (i, e) in ind.genome.iter().enumerate() {
                r_genome[*e] = i;
            }

            let mut valid = 0;
            for rule in &self.rules {
                let a = rule.0;
                let b = rule.1;
                let c = rule.2;
                if      ((r_genome[a] < r_genome[b]) && (r_genome[b] < r_genome[c]))
                    ||  ((r_genome[a] > r_genome[b]) && (r_genome[b] > r_genome[c]))
                    {
                    valid += 1;
                }
            }
            valid as f64 / self.rules.len() as f64
        }
    }

fn solve_solvable(elements: usize, triples: usize, generations: usize, unchanged: usize, popsize: usize, mutation: &str) {
    // Permutation erzeugen
    let solution = Individual::new_random(elements);
    // Regeln erzeugen
    let problem = BTWNSProblem::new_solvable_from_individual(&solution, elements, triples);
    println!("Lösung: {}", solution);
    println!("Regeln: {}", problem);

    // Evoli
    let mut individuals = Vec::new();
    for _ in 0..popsize {
        individuals.push(Individual::new_random(elements));
    }
    let mut ind = Individual::new_random(elements);
    let mut max = 0.0;
    for i in 0..generations {
        let rating = problem.rate(&ind);
        if rating == 1.0 {
            println!("{}", "=".repeat(80));
            println!("gefundene Lösung ({} Durchläufe): {}", i, ind);
            max = rating;
            break;
        }
        ind.mutate();
        max = max.max(rating);

    }
    println!("Max: {}", max);

}
