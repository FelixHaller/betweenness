extern crate rand;
extern crate clap;
extern crate roulette;
extern crate random_choice;

use clap::{Arg, App, SubCommand};
use std::fmt;
use rand::Rng;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::prelude::*;
use std::error::Error;
use std::path::Path;

use roulette::Roulette;
use self::random_choice::random_choice;

const DEFAULT_MUTATION: &str = "swap";
const DEFAULT_ELEMENTS: &str = "100";
const DEFAULT_TRIPLES: &str = "50";
const DEFAULT_GENERATIONS: &str = "1000";
const DEFAULT_MAXUNCHANGED: &str = "500";
const DEFAULT_POPULATIONSIZE: &str = "1000";

//const PARENTSELFAC: f64 = 0.5;
const MUTATIONRATE: f64 = 1.0;

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
        .about("Dieses Programm bearbeitet das Problem MAXIMUM BETWEENNESS mittels evolutionärem Algorithmus.")
        .author("Felix Haller <felix.haller@stud.htwk-leipzig.de>")
        .arg(Arg::with_name("verbose")
            .short("v")
            .long("verbose")
            .help("Mehr Ausgaben."))
        .arg(Arg::with_name("mutation")
            .short("m")
            .value_name("MUTATION")
            .help("Mutationsstrategie")
            .possible_values(&["swap", "shift", "invert"])
            .takes_value(true)
            .required(true)
            .default_value(DEFAULT_MUTATION))
        .arg(Arg::with_name("elements")
            .short("e")
            .value_name("ELEMENTS")
            .help("Anzahl der Elemente der Permutation")
            .takes_value(true)
            .required(true)
            .validator(is_usize)
            .default_value(DEFAULT_ELEMENTS))
        .arg(Arg::with_name("triples")
            .short("t")
            .value_name("TRIPLES")
            .help("Anzahl der Triples")
            .takes_value(true)
            .required(true)
            .validator(is_usize)
            .default_value(DEFAULT_TRIPLES))
        .arg(Arg::with_name("generations")
            .short("g")
            .value_name("GENERATIONS")
            .help("Anzahl max. Generationen")
            .takes_value(true)
            .required(true)
            .validator(is_usize)
            .default_value(DEFAULT_GENERATIONS))
        .arg(Arg::with_name("maxunchanged")
            .short("u")
            .value_name("MAXUNCHANGED")
            .help("Anzahl der max. aufeinanderfolgenden Generationen ohne Verbesserung")
            .takes_value(true)
            .required(true)
            .validator(is_usize)
            .default_value(DEFAULT_MAXUNCHANGED))
        .arg(Arg::with_name("popsize")
            .short("p")
            .value_name("POPULATIONSIZE")
            .help("Populationsgröße festlegen")
            .takes_value(true)
            .required(true)
            .validator(is_even_usize)
            .default_value(DEFAULT_POPULATIONSIZE))
        .arg(Arg::with_name("filename")
            .short("f")
            .value_name("FILENAME")
            .help("Triple aus Datei laden")
            .takes_value(true)
            .required(false)
            .conflicts_with("elements")
            .conflicts_with("triples"))
        .subcommand(SubCommand::with_name("hillclimber")
            .about("Versuch der Lösung eines Problems mittels Hillclimber.")
            .version("1.0")
            .author("Felix Haller <felix.haller@stud.htwk-leipzig.de>")
            .arg(Arg::with_name("verbose")
                .short("v")
                .help("Mehr Ausgaben."))
            .arg(Arg::with_name("filename")
                .short("f")
                .value_name("FILENAME")
                .help("Triple aus Datei laden")
                .takes_value(true)
                .required(true))
            .arg(Arg::with_name("iterations")
                .short("i")
                .value_name("ITERATIONS")
                .help("Anzahl max. Iterationen.")
                .takes_value(true)
                .required(true)
                .validator(is_usize)
                .default_value(DEFAULT_GENERATIONS))
        )
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
    let max_unchanged = matches.value_of("maxunchanged").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let popsize = matches.value_of("popsize").unwrap().parse::<usize>().expect("Fehler beim Parsen");
    let verbose = matches.is_present("verbose");

    if let Some(matches) = matches.subcommand_matches("hillclimber") {
        let verbose = matches.is_present("verbose");
        let iterations = matches.value_of("iterations").unwrap().parse::<usize>().expect("Fehler beim Parsen");
        if let Some(filename) = matches.value_of("filename") {
            solve_with_hillclimber_from_file(filename, iterations,verbose);
        }
    }

    if let Some(filename) = matches.value_of("filename") {
        solve_from_file(filename, generations, max_unchanged, popsize, mutation, verbose);
    }
    else {
        solve_solvable(elements, triples, generations, max_unchanged, popsize, mutation, verbose);
    }




}
#[derive(Debug)]
struct Individual {
    //Generation in der das Individuum erschaffen wurde
    generation: usize,
    // Generation der Mutter
    mother_gen: usize,
    //Generation des Vaters
    father_gen: usize,
    // Das Genome
    genome: Vec<usize>,
    // Die Bewertung
    rating: f64,
}

    impl fmt::Display for Individual {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{} | Rating: [{}] | Generation:{}, Elterngenerationen:{}, {}", format!("{:?}", self.genome), self.rating, self.generation,  self.mother_gen, self.father_gen)
        }
    }

    impl Individual {
        fn new_dummy() -> Individual {
            Individual{ generation: 0, mother_gen: 0, father_gen: 0, genome: Vec::new(), rating: 0.0}
        }
        fn new_random(elements: usize) -> Individual {
            let mut genome = (0..elements).collect::<Vec<usize>>();
            rand::thread_rng().shuffle(&mut genome);
            Individual{ generation: 0, mother_gen: 0, father_gen: 0, genome: genome, rating: 0.0}
        }
        fn new_random_from_elements(elements: &BTreeSet<usize>) -> Individual {
            let mut genome = elements.iter().cloned().collect::<Vec<usize>>();
            rand::thread_rng().shuffle(&mut genome);
            Individual{ generation: 0, mother_gen: 0, father_gen: 0, genome: genome, rating: 0.0}
        }

        fn get_copy(&self) -> Individual {
            Individual{ generation: self.generation, mother_gen: self.mother_gen, father_gen: self.father_gen, genome: self.genome.to_vec(), rating: self.rating}
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
            let a = rng.gen_range(0, self.genome.len()-1);
            let b = rng.gen_range(a+1, self.genome.len());

            self.genome[a..b].reverse();

        }

        /// verschiebende Mutation
        fn mutate_shift(&mut self) {
            //(gut geeignet laut Skript)
            let mut rng = rand::thread_rng();
            let len = self.genome.len();
            let e = self.genome.remove(rng.gen_range(0, len));
            self.genome.insert(rng.gen_range(0, len), e);

        }

        fn recombine(&self, other : &Individual, generation: usize) -> (Individual, Individual) {
            // ORDNUNGSREKOMBINATION (ohne Klone) (Buch S.29)
            let crosspoint = rand::thread_rng().gen_range(1, self.genome.len()-1);

            let mut child1 = Individual{ generation, mother_gen: self.generation, father_gen: other.generation, genome: self.genome[0..crosspoint].to_vec(), rating: 0.0};
            let mut child2 = Individual{ generation, mother_gen: self.generation, father_gen: other.generation, genome: other.genome[0..crosspoint].to_vec(), rating: 0.0};

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
        fn new_solvable_from_individual(solution: &Individual, triples: usize) -> BTWNSProblem {
            let mut rules = BTreeSet::<(usize,usize,usize)>::new();
            let mut rng = rand::thread_rng();
            while rules.len() < triples {
                let i_c = rng.gen_range(2, solution.genome.len());
                let i_b = rng.gen_range(1, i_c);
                let i_a = rng.gen_range(0, i_b);
                let c = solution.genome[i_c];
                let b = solution.genome[i_b];
                let a = solution.genome[i_a];
                rules.insert((a,b,c));
            }
            BTWNSProblem{rules : rules}
        }

        fn new_from_file(filename: &str) -> BTWNSProblem {
            let mut rules = BTreeSet::new();
            let path = Path::new(filename);
            let mut f = match File::open(filename) {
                // The `description` method of `io::Error` returns a string that
                // describes the error
                Err(why) => panic!("Konnte Datei {} nicht öffnen: {}", path.display(),
                                   why.description()),
                Ok(file) => file,
            };

            let mut contents = String::new();
            f.read_to_string(&mut contents).expect("something went wrong reading the file");
            //println!("{}", contents);

            let lines: Vec<&str> = contents.split('\n').collect();
            let rule_line = lines[0].replace(" ", "").replace("{", "").replace("}", "");
            let tuples = rule_line.split("),").collect::<Vec<&str>>();
            for tuple in tuples{
                let clean = tuple.replace("(", "").replace(")", "");
                let abc: Vec<&str> = clean.split(",").collect();
                let a = abc[0].parse::<usize>().expect("Fehler in der Eingabedatei");
                let b = abc[1].parse::<usize>().expect("Fehler in der Eingabedatei");
                let c = abc[2].parse::<usize>().expect("Fehler in der Eingabedatei");
                rules.insert((a,b,c));
            }
            BTWNSProblem{rules}
        }

        fn rate(&self, ind: &Individual, print_invalid: bool) -> f64 {
            // genome und r_genome erzeugen um später in konstanter Zeit die Position eines
            // Elements zu erfahren
            let mut r_genome = vec![0; ind.genome.iter().cloned().max().expect("Individuum ohne Genom")+1];
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
                else if print_invalid {
                    eprint!("({}, {}, {}) ", a, b, c);
                }
            }
            valid as f64 / self.rules.len() as f64
        }
    }

fn _solve(problem: BTWNSProblem, elements: BTreeSet<usize>, generations: usize, max_unchanged: usize, popsize: usize, mutation: MutationStrategy, verbose: bool){
    println!("Regeln: {}", problem);
    println!("Elemente (effektiv): {}", elements.len());
    let mut best_individual = Individual::new_dummy();
    let mut rng = rand::thread_rng();

    //Start-Population zufällig erzeugen
    let mut parent_population = Vec::new();
    for _ in 0..popsize {
        parent_population.push(Individual::new_random_from_elements(&elements));
    }

    // Initiales Bewerten
    for individual in parent_population.iter_mut() {
        //Terminierungsbedingung prüfen
        individual.rating = problem.rate(&individual, false);
        if individual.rating == 1.0 {
            exit(format!("Lösung (in Startpopulation): {}", individual).as_str(), 0);
        }
        if individual.rating > best_individual.rating {
            best_individual = individual.get_copy();
        }


    }
    //Terminierungsbedingung 1 (Anzahl Generationen)
    for generation in 0..generations {
        let avg_rating = parent_population.iter().fold(0.0, |mut sum, val| {sum += val.rating; sum}) / popsize as f64;
        let old_best_rating = best_individual.rating;
        let mut child_population = Vec::new();

        //Terminierungsbedingung 2 prüfen (inaktiv)
        if false {
            if generation - best_individual.generation > max_unchanged {
                exit(format!("Keine Verbesserung seit {} Generationen.\nBestes Individuum bisher: {}", max_unchanged, best_individual).as_str(), 1);
            }
        }

        //Mischen
        rng.shuffle(&mut parent_population);


        // Elternselektion
        // Alle Individuen sollten Chancen erhalten
        // Fitness Proportionale Selektion oder Stochastisches Universelles Sampling
//        let mut selected_parents = stoch_univ_samp(&parent_population, popsize/2);
        let mut selected_parents = fit_prop_sel(&parent_population, popsize/2, &mut rng);

        // Rekombination
        while let Some(individual) = selected_parents.pop() {
            let partner = selected_parents.pop().expect("Populationszahl war ungerade!");
            let (child1, child2) = individual.recombine(&partner, generation);
            let (child3, child4) = individual.recombine(&partner, generation);
            child_population.push(child1);
            child_population.push(child2);
            child_population.push(child3);
            child_population.push(child4);
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
            individual.rating = problem.rate(&individual, false);
            // Terminierungsbedingung "Optimum" prüfen
            if individual.rating == 1.0 {
                exit(format!("Lösung (nach Generation {}):\n{}", generation, individual).as_str(), 0);
            }
            if individual.rating > best_individual.rating {
                best_individual = individual.get_copy();
            }
        }
        if verbose {
            if best_individual.rating > old_best_rating || generation == 0 {
                println!("[{}] {}", generation, best_individual);
            }
        }
        else {
            if best_individual.rating > old_best_rating  || generation == 0 {
                println!("[{}]: {:.*} (ø{:.*})", generation, 4, best_individual.rating, 4, avg_rating);
//                println!("({}, {:.*} {:.*})", generation, 4, best_individual.rating, 4, avg_rating);
                println!("---");
            }

        }
        eprint!("[{}]: {:.*} (ø{:.*})\r", generation, 4, best_individual.rating, 4, avg_rating);
//        eprint!("[{}]: {:.*}\r", generation, 4, best_individual.rating);

        // Umweltselektion
        parent_population.append(&mut child_population);

        // Q-Stufige zweifache Turnierselektion? (Eltern gegen Kinder oder Gegener jeweils innerhalb der Gruppen?)
        // ("alte" Individuen, verschwinden meist automatisch mit der Zeit, wenn nach Wahrscheinlichkeiten aussortiert wird)
        parent_population = q_step_tournament_sel(parent_population, popsize, 3, &mut rng);

        // rangbasierte Selektion?

    }
    eprint!("Nicht erfüllte Regeln: ");
    problem.rate(&best_individual, true);
    eprintln!("");
    eprintln!("Keine Lösung gefunden in {} Generationen.", generations);
    eprintln!("Bestes Individuum bisher: {}", best_individual);
}

fn solve_from_file(filename: &str, generations: usize, max_unchanged: usize, popsize: usize, mutation: MutationStrategy, verbose: bool) {
    let problem = BTWNSProblem::new_from_file(filename);
    println!("{} Triple (Regeln) geladen", problem.rules.len());
    let mut elements = BTreeSet::new();

    for triple in problem.rules.iter() {
        elements.insert(triple.0);
        elements.insert(triple.1);
        elements.insert(triple.2);
    }
    println!("Enthaltene Elemente: {:?} (Anzahl: {})", elements, elements.len());
    _solve(problem, elements, generations, max_unchanged, popsize, mutation, verbose);

}

fn solve_solvable(elements_n: usize, triples_m: usize, generations: usize, max_unchanged: usize, popsize: usize, mutation: MutationStrategy, verbose: bool) {
    // Permutation erzeugen
    let solution = Individual::new_random(elements_n);
    println!("mögliche Lösung: {:?}", solution.genome);
    // Regeln erzeugen
    let problem = BTWNSProblem::new_solvable_from_individual(&solution, triples_m);

    let mut elements = BTreeSet::new();
    for triple in problem.rules.iter() {
        elements.insert(triple.0);
        elements.insert(triple.1);
        elements.insert(triple.2);
    }

    _solve(problem, elements, generations, max_unchanged, popsize, mutation, verbose);

}

fn solve_with_hillclimber_from_file(filename: &str, runs: usize, verbose: bool){
    let problem = BTWNSProblem::new_from_file(filename);
    println!("{} Triple (Regeln) geladen", problem.rules.len());
    let mut elements = BTreeSet::new();

    for triple in problem.rules.iter() {
        elements.insert(triple.0);
        elements.insert(triple.1);
        elements.insert(triple.2);
    }
    println!("Enthaltene Elemente: {:?}", elements);
    println!("Elemente effektiv: {}", elements.len());
    let mut climber = Individual::new_random_from_elements(&elements);
    for i in 0..runs {
        let rating = problem.rate(&climber, false);
        climber.rating = rating;
        if rating == 1.0 {
            exit(format!("Lösung (nach Durchlauf {}):\n{}", i, climber).as_str(), 0);
        }
        let mut new_climber = Individual::new_random_from_elements(&elements);
        new_climber.generation = i;
        new_climber.father_gen = i-1;
        new_climber.mother_gen = i-1;
        let new_rating = problem.rate(&new_climber, false);
        new_climber.rating = new_rating;
        if new_climber.rating > rating {
            climber = new_climber;
        }
        eprint!("[{}]: {:.*}\r", i, 4, climber.rating);
        if new_rating > rating || i == 0 {
            if verbose {
                println!("[{}] | {:.*}", i, 4, climber);
                println!("---");
            }
            else {
                println!("[{}]: {:.*}", i, 4, climber.rating);
            }
        }

    }
    eprintln!("Keine Lösung gefunden in {} Durchläufen.", runs);
    eprintln!("Bestes Individuum bisher: {}", climber);
    exit("Ende", 0);

}

fn q_step_tournament_sel(population : Vec<Individual>, to_select: usize, q: u8, rng: &mut rand::ThreadRng) -> Vec<Individual>{
    let mut results = Vec::new();
    let mut champions = Vec::new();
    for individual in population.iter() {
        let mut wins = 0;
        for _ in 0..q {
            if individual.rating > rng.choose(&population).expect("Population war leer!").rating {
                wins += 1;
            }
        }
        results.push((
            individual.get_copy(),
            wins));
    }
    results.sort_by(|a, b| a.1.cmp(&b.1).reverse());
    for i in 0..to_select {
        let result = &results[i];
        champions.push(
            result.0.get_copy());
    }
    champions
}

fn stoch_univ_samp(parent_population: &Vec<Individual>, to_select: usize) -> Vec<Individual> {

    let weights: Vec<f64> = parent_population.iter().map(|x| x.rating).collect();
    let choices = random_choice().random_choice_f64(&parent_population, &weights, to_select);
    let mut selected_parents = Vec::new();

    for choice in choices {
        selected_parents.push(choice.get_copy());
    }

    selected_parents
}

fn fit_prop_sel(parent_population: &Vec<Individual>, to_select: usize, mut rng: &mut rand::ThreadRng) -> Vec<Individual> {
    let mut parents = Vec::new();
    let mut selected_parents = Vec::new();

    for individual in parent_population.iter() {
        let rating = individual.rating;
        parents.push((individual.get_copy(),rating));
    }
    let roulette = Roulette::new(parents);
    for _ in 0..to_select {
        let chosen = roulette.next(&mut rng);
        selected_parents.push(chosen.get_copy());
    }
    selected_parents

}
