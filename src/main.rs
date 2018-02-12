extern crate rand;
extern crate clap;
extern crate roulette;

use clap::{Arg, App};
use std::fmt;
use rand::Rng;
use std::collections::BTreeSet;
use std::collections::HashMap;

use roulette::Roulette;

const DEFAULT_MUTATION: &str = "swap";
const DEFAULT_ELEMENTS: &str = "100";
const DEFAULT_TRIPLES: &str = "50";
const DEFAULT_GENERATIONS: &str = "1000";
const DEFAULT_MAXUNCHANGED: &str = "500";
const DEFAULT_POPULATIONSIZE: &str = "1000";


const PARENTSELFAC: f64 = 0.5;
const MUTATIONRATE: f64 = 0.3;

enum MutationStrategy {
    Swap,
    Shift,
    Invert
}

fn exit(msg: &str, status: i32) {
    eprintln!("{}", msg);
    std::process::exit(status);
}

fn is_usize(s: String) -> Result<(), String> {
    match s.parse::<usize>() {
        Ok(_) => Ok(()),
        Err(_) => Err("Ist kein Integer-Wert".to_string())
    }
}

fn is_even_usize(s: String) -> Result<(), String> {
    match s.parse::<usize>() {
        Ok(u) =>  match u%2 {
            0 => Ok(()),
            _ => Err("Dieser Algorithmus akzeptiert nur gerade Werte".to_string())
        },
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
            .default_value(DEFAULT_MUTATION)
            .takes_value(true)
            .required(false))
        .arg(Arg::with_name("elements")
            .short("e")
            .value_name("ELEMENTS")
            .help("Anzahl der Elemente der Permutation")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value(DEFAULT_ELEMENTS))
        .arg(Arg::with_name("triples")
            .short("t")
            .value_name("TRIPLES")
            .help("Anzahl der Triples")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value(DEFAULT_TRIPLES))
        .arg(Arg::with_name("generations")
            .short("g")
            .value_name("GENERATIONS")
            .help("Anzahl max. Generationen")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value(DEFAULT_GENERATIONS))
        .arg(Arg::with_name("maxunchanged")
            .short("u")
            .value_name("MAXUNCHANGED")
            .help("Anzahl der max. aufeinanderfolgenden Generationen ohne Verbesserung")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value(DEFAULT_MAXUNCHANGED))
        .arg(Arg::with_name("popsize")
            .short("p")
            .value_name("POPULATIONSIZE")
            .help("Populationsgröße festlegen")
            .takes_value(true)
            .required(false)
            .validator(is_even_usize)
            .default_value(DEFAULT_POPULATIONSIZE))
        .get_matches();

    let mutation = match matches.value_of("mutation").expect("Fehler beim Parsen"){
        "swap" => MutationStrategy::Swap,
        "shift" => MutationStrategy::Shift,
        "invert" => MutationStrategy::Invert,
        _ => panic!("Fehler beim Parsen")
    };
    let elements = matches.value_of("elements").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let triples = matches.value_of("triples").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let generations = matches.value_of("generations").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let maxunchanged = matches.value_of("maxunchanged").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let popsize = matches.value_of("popsize").unwrap().parse::<usize>().expect("Fehler beim Parsen");

    solve_solvable(elements, triples, generations, maxunchanged, popsize, mutation);
}
#[derive(Debug)]
struct Individual {
    generation: usize,
    // Das Genome
    genome: Vec<usize>,
    // Die Bewertung
    rating: f64,
}

    impl fmt::Display for Individual {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[{}]: {}", self.generation, format!("{:?}", self.genome))
        }
    }

    impl Individual {
        fn new_random(elements: usize) -> Individual {
            let mut genome = (0..elements).collect::<Vec<usize>>();
            rand::thread_rng().shuffle(&mut genome);
            Individual{generation : 1, genome: genome, rating: 0.0}
        }

        fn mutate_swap(&mut self) {
            // vertauschende Mutation
            let mut rng = rand::thread_rng();
            let a = rng.gen_range(0, self.genome.len());
            let b = rng.gen_range(0, self.genome.len());

            self.genome.swap(a, b);
        }

        /// invertierende Mutation
        fn mutate_invert(&mut self) {
            let mut rng = rand::thread_rng();
            let a = rng.gen_range(0, self.genome.len());
            let b = rng.gen_range(0, self.genome.len());

            self.genome[a..(b+1)].reverse();
        }

        /// verschiebende Mutation
        fn mutate_shift(&mut self) {
            //(gut geeignet laut Skript)
            let mut rng = rand::thread_rng();
            let len = self.genome.len();
            let e = self.genome.remove(rng.gen_range(0, len));
            self.genome.insert(rng.gen_range(0, len), e);

        }

        fn recombine(&self, other : Individual) -> (Individual, Individual) {
            // ORDNUNGSREKOMBINATION (ohne Klone) (Buch S.29)
            let new_gen = self.generation + 1;
            let crosspoint = rand::thread_rng().gen_range(1, self.genome.len()-1);

            let mut child1 = Individual{generation: new_gen, genome: self.genome[0..crosspoint].to_vec(), rating: 0.0};
            let mut child2 = Individual{generation: new_gen, genome: other.genome[0..crosspoint].to_vec(), rating: 0.0};

            for parent_genome in vec![&self.genome, &other.genome] {
                let mut genome = parent_genome.to_vec(); //TODO: ohne Kopie lösen?
                while let Some(element) = genome.pop() {
                    if ! child1.genome.contains(&element) {
                        child1.genome.push(element);
                    }
                    if ! child2.genome.contains(&element) {
                        child2.genome.push(element);
                    }
                }
            }
            (child1, child2)
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
            let mut rng = rand::thread_rng();
            while rules.len() < triples {
                let i_c = rng.gen_range(2, elements);
                let i_b = rng.gen_range(1, i_c);
                let i_a = rng.gen_range(0, i_b);
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

fn solve_solvable(elements: usize, triples: usize, generations: usize, max_unchanged: usize, popsize: usize, mutation: MutationStrategy) {
    // Permutation erzeugen
    let solution = Individual::new_random(elements);
    // Regeln erzeugen
    let problem = BTWNSProblem::new_solvable_from_individual(&solution, elements, triples);
    println!("Lösung: {}", solution);
    println!("Regeln: {}", problem);

    //+++++++++++++++++++++++++++++++++
    //Zähler für wiederholte beste Werte (Maximum, inGeneration)
    struct MaxQuality {
        value: f64,
        generation: usize,
    }
    let mut max_quality = MaxQuality{value: 0.0, generation: 0};
    let mut rng = rand::thread_rng();

    //Start-Population zufällig erzeugen
    let mut parent_population = Vec::new();
    for _ in 0..popsize {
        parent_population.push(Individual::new_random(elements));
    }

    // Initiales Bewerten
    for individual in parent_population.iter_mut() {
        //Terminierungsbedingung prüfen
        individual.rating = problem.rate(&individual);
        if individual.rating == 1.0 {
            exit(format!("Lösung: {}", individual).as_str(), 0);
        }
        if individual.rating > max_quality.value {
            max_quality = MaxQuality{value: individual.rating, generation: 0};
        }


    }
    //Terminierungsbedingung 1 (Anzahl Generationen)
    for generation in 1..generations {
        let mut child_population = Vec::new();

        // Terminierungsbedingung 2 prüfen
//        if generation - max_quality.generation > max_unchanged {
//            exit(format!("Keine Verbesserung seit {} Generationen", max_unchanged).as_str(), 1);
//        }


        rng.shuffle(&mut parent_population);
        

        // Elternselektion
        // Alle Individuen sollten Chancen erhalten
        // Fitness Proportionale Selektion
        let to_select = popsize/2;
        let mut parents = Vec::new();
        let mut selected_parents = Vec::new();

        for individual in parent_population.iter() {
            let rating = individual.rating;
            parents.push((
                Individual {
                generation: individual.generation,
                genome: individual.genome.to_vec(),
                rating: individual.rating,
            }, rating));
        }
        let roulette = Roulette::new(parents);
        for _ in 0..to_select {
            let chosen = roulette.next(&mut rng);
//            println!("chosen: {}", chosen);
            selected_parents.push(
                Individual {
                    generation: chosen.generation,
                    genome: chosen.genome.to_vec(),
                    rating: chosen.rating,
                }
            );
        }


        // Rekombination
        while let Some(individual) = selected_parents.pop() {
            let (child1, child2) = individual.recombine(selected_parents.pop().expect("Populationszahl war ungerade!"));
            child_population.push(child1);
            child_population.push(child2);
        }


        // Mutation der Kinder
        for individual in child_population.iter_mut() {
            if rng.gen_range(0.0, 1.0) < MUTATIONRATE {
                match mutation {
                    MutationStrategy::Swap => individual.mutate_swap(),
                    MutationStrategy::Shift => individual.mutate_shift(),
                    MutationStrategy::Invert => individual.mutate_invert(),

                }
            }
        }

        // Bewertung und auf Optimum prüfen
        for individual in child_population.iter_mut() {
            individual.rating = problem.rate(&individual);
            // Terminierungsbedingung "Optimum" prüfen
            if individual.rating == 1.0 {
                exit(format!("Lösung: {}", individual).as_str(), 0);
            }
            if individual.rating > max_quality.value {
                max_quality = MaxQuality{value: individual.rating, generation};
            }
        }
        println!("Generation: {} \t\t max: {}", generation, max_quality.value);


        // Umweltselektion
        //Eltern werden durch Kinder vollständig ersetzt
        parent_population = child_population;
        // gleiche Individuen aussortieren
//        parent_population.sort_by(|a, b| a.genome.cmp(&b.genome));
//        parent_population.dedup_by(|a, b| a.genome.eq(&b.genome));


        // evtl. Q-Stufige Turnierselektion?
        // rangbasierte Selektion
        // ("alte" Individuen, verschwinden meist automatisch mit der Zeit, wenn nach Wahrscheinlichkeiten aussortiert wird)

    }
}

