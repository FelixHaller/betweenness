extern crate rand;
extern crate clap;

use clap::{Arg, App};
use std::fmt;
use rand::Rng;
use std::collections::BTreeSet;

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
        .arg(Arg::with_name("maxunchanged")
            .short("u")
            .value_name("MAXUNCHANGED")
            .help("Anzahl der max. aufeinanderfolgenden Generationen ohne Verbesserung")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .default_value("500"))
        .arg(Arg::with_name("popsize")
            .short("p")
            .long("pop")
            .value_name("POPULATIONSIZE")
            .help("Populationsgröße festlegen")
            .takes_value(true)
            .required(false)
            .validator(is_usize)
            .validator(is_even_usize)
            .default_value("100"))
        .get_matches();

    let mutation = matches.value_of("mutation").expect("Fehler beim Parsen");
    let elements = matches.value_of("elements").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let triples = matches.value_of("triples").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let generations = matches.value_of("generations").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let maxunchanged = matches.value_of("maxunchanged").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let popsize = matches.value_of("popsize").unwrap().parse::<usize>().expect("Fehler beim Parsen");

    solve_solvable(elements, triples, generations, maxunchanged, popsize, mutation);
}

struct Individual {
    generation: usize,
    // Das Genome
    genome: Vec<usize>,
}

    impl fmt::Display for Individual {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "[{}]: {}", self.generation, format!("{:?}", self.genome))
        }
    }

    impl Individual {
        fn new_random(elements: usize) -> Individual {
            let mut genome = Vec::new();
            for i in 0..elements {
                genome.push(i);
            }
            rand::thread_rng().shuffle(&mut genome);


            Individual{generation : 1, genome: genome}
        }

        fn mutate(&mut self, mutation: &str) {

            // vertauschende Mutation
            let a = rand::thread_rng().gen_range(0, self.genome.len());
            let b = rand::thread_rng().gen_range(0, self.genome.len());

            self.genome.swap(a, b);


            // invertierende Mutation

            // verschiebende Mutation (gut geeignet laut Skript)

        }
        fn recombine(&self, other : Individual) -> (Individual, Individual) {
            // ORDNUNGSREKOMBINATION (ohne Klone) (Buch S.29)
            let new_gen = self.generation + 1;
            let crosspoint = rand::thread_rng().gen_range(1, self.genome.len()-1);

            let mut child1 = Individual{generation: new_gen, genome: self.genome[0..crosspoint].to_vec()};
            let mut child2 = Individual{generation: new_gen, genome: other.genome[0..crosspoint].to_vec()};

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

fn solve_solvable(elements: usize, triples: usize, generations: usize, max_unchanged: usize, popsize: usize, mutation: &str) {
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


    //Start-Population zufällig erzeugen
    let mut population = Vec::new();
    for _ in 0..popsize {
        population.push(Individual::new_random(elements));
    }

    // Initiales Bewerten
    let mut ratings = Vec::new();
    for individual in population.iter() {
        //Terminierungsbedingung prüfen
        let rating = problem.rate(&individual);
        if rating == 1.0 {
            exit(format!("Lösung: {}", individual).as_str(), 0);
        }
        if rating > max_quality.value {
            max_quality = MaxQuality{value: rating, generation: 0};
        }
        ratings.push(rating);

    }
    for generation in 1..generations {
        // Terminierungsbedingungen prüfen
//        if generation - max_quality.generation > max_unchanged {
//            exit(format!("Keine Verbesserung seit {} Generationen", max_unchanged).as_str(), 1);
//        }
        if max_quality.value == 1.0 {
            for individual in population.iter() {
                if problem.rate(&individual) == 1.0 {
                    exit(format!("Lösung: {}", individual).as_str(), 0);
                }
            }
        }
        // Paarungsselektion
        //SKIP

        // Rekombination
        let mut new_population = Vec::new();
        rand::thread_rng().shuffle(&mut population);
        while let Some(individual) = population.pop() {
            let (child1, child2) = individual.recombine(population.pop().expect("Populationszahl war ungerade!"));
            new_population.push(child1);
            new_population.push(child2);
        }
        //Eltern werden durch Kinder vollständig ersetzt
        population = new_population;

        // Mutation der Kinder
        for individual in population.iter_mut() {
            individual.mutate(mutation);
        }


        // Bewertung
        ratings = Vec::new();
        for individual in population.iter() {
            let rating = problem.rate(&individual);
            if rating > max_quality.value {
                max_quality = MaxQuality{value: rating, generation};
            }
            ratings.push(rating);
        }
        println!("Generation: {} \t\t max: {}", generation, max_quality.value);

        // Umweltselektion
    }


//OLD STUFF


//    let mut max = 0.0;
//    for i in 0..generations {
//        let rating = problem.rate(&ind);
//        if rating == 1.0 {
//            println!("{}", "=".repeat(80));
//            println!("gefundene Lösung ({} Durchläufe): {}", i, ind);
//            max = rating;
//            break;
//        }
//        ind.mutate();
//        max = max.max(rating);
//
//    }
//    println!("Max: {}", max);

}
