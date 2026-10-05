use crate::alphabet::alphabet::Alphabet;
use crate::seq::mutation::MutationEvent;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MutationClass {
  Substitution,
  Ambiguous,
  Indel,
}

pub fn classify_mutation(event: &MutationEvent, alphabet: &Alphabet) -> MutationClass {
  match event {
    MutationEvent::Substitution(sub) if alphabet.is_canonical(sub.reff()) && alphabet.is_canonical(sub.qry()) => {
      MutationClass::Substitution
    },
    MutationEvent::Substitution(_) => MutationClass::Ambiguous,
    MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => MutationClass::Indel,
  }
}
