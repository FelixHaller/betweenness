extern crate rand;

use std::fmt;
use rand::Rng;
use std::collections::HashSet;

/// Anzahl der Zeichen
const ELEMENTS: usize = 10;
/// Anzahl der Contraints
const TRIPLES: usize = 50;
/// max. Anzahl der Versuche
const MAXGENERATIONS: usize = 100000;
/// max. Anzahl der aufeinanderfolgenden Genrationen ohne Verbesserung
const MAXUNCHANGED: usize = 100;
/// Individuen in Population
const POPULATION: usize = 100;

type t_genome = [usize; ELEMENTS];

struct Individual {
    id : u32,
    generation: u32,
    // Das Genome
    genome: t_genome,
    // Eine umgekehrte Zuordnung um schnell den Index eines Elements zu erfahren
    r_genome: t_genome,
}

    impl fmt::Display for Individual {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "#{}[{}]: {}", self.id, self.generation, format!("{:?}", self.genome))
        }
    }

    impl Individual {
        fn new_random() -> Individual {
            let mut genome = [0usize; ELEMENTS];
            for i in 0..ELEMENTS {
                genome[i] = i;
            }
            rand::thread_rng().shuffle(&mut genome);


            Individual{id : 0, generation : 1, genome: genome, r_genome: Individual::create_rgenome(genome)}
        }

        fn create_rgenome(genome : t_genome) -> t_genome {
            let mut r_genome = [0usize; ELEMENTS];
            for (i, element) in genome.iter().enumerate() {
                r_genome[*element] = i;
            }
            r_genome
        }

        fn refresh_rgenome(&mut self) {
            self.r_genome = Individual::create_rgenome(self.genome);
        }

        fn pos2el(&self, pos: usize) -> usize {
            self.genome[pos]
        }

        fn el2pos(&self, element: usize) -> usize {
            self.r_genome[element]
        }

        fn mutate(&mut self) {
            println!("Vorher: {}", self);

            // vertauschende Mutation
            let a = rand::thread_rng().gen_range(0, self.genome.len());
            let b = rand::thread_rng().gen_range(0, self.genome.len());

            let s = self.genome[a];
            self.genome[a] = self.genome[b];
            self.genome[b] = s;

            self.refresh_rgenome();

            println!("Nachher: {}", self);
            // invertierende Mutation

            // verschiebende Mutation (gut geeignet laut Skript)

        }
        fn recombine(&self, other : Individual) {
            // ORDNUNGSREKOMBINATION (Buch S.29)
        }
    }

struct BTWNSProblem {
    rules : HashSet<(usize,usize,usize)>,
}

    impl fmt::Display for BTWNSProblem {
        fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "{:?}", self.rules)
        }
    }

    impl BTWNSProblem {
        fn new_solvable_from_individual(solution: &Individual) -> BTWNSProblem {
            let mut rules : HashSet<(usize,usize,usize)> = HashSet::new();
            while rules.len() < TRIPLES {
                let i_c = rand::thread_rng().gen_range(2, ELEMENTS);
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
            let mut valid = 0;
            for rule in &self.rules {
                let a = rule.0;
                let b = rule.1;
                let c = rule.2;
                if      ((ind.el2pos(a) < ind.el2pos(b)) && (ind.el2pos(b) < ind.el2pos(c)))
                    ||  ((ind.el2pos(a) > ind.el2pos(b)) && (ind.el2pos(b) > ind.el2pos(c)))
                    {
                    valid += 1;
                }
            }
            valid as f64 / self.rules.len() as f64
        }
    }

fn solve_solvable() {
    // Permutation erzeugen
    let solution = Individual::new_random();
    // Regeln erzeugen
    let problem = BTWNSProblem::new_solvable_from_individual(&solution);
    println!("Lösung: {:?}", solution.genome);
    println!("Regeln: {}", problem);

    // Ab hier Evoli
    let mut ind = Individual::new_random();
    let mut max = 0.0;
    for x in 0..MAXGENERATIONS {
        let rating = problem.rate(&ind);
        if rating == 1.0 {
            println!("========================================");
            println!("Regeln: {}", problem);
            println!("Lösung: {}", solution);
            println!("gefundene Lösung: {}", ind);
            break;
        }
        ind.mutate();
        println!("{}: {}", x, rating);
        if rating > max {
            max = rating;
        }

    }
    println!("Max: {}", max);

}



//fn solve(problem : BTWNSProblem) {
//
//}



fn main() {
    solve_solvable();
}
