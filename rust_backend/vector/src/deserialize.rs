use serde::de::{Deserializer, Error, SeqAccess, Visitor};
use serde::Deserialize;
use std::fmt;
use std::marker::PhantomData;
use {Field, LorentzVector};

struct LorentzVectorVisitor<T: Field> {
    _marker: PhantomData<fn() -> LorentzVector<T>>,
}

impl<'de, T: Field + Deserialize<'de>> Visitor<'de> for LorentzVectorVisitor<T> {
    type Value = LorentzVector<T>;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("four floats")
    }

    fn visit_seq<M>(self, mut access: M) -> Result<Self::Value, M::Error>
    where
        M: SeqAccess<'de>,
    {
        let t = access
            .next_element::<T>()?
            .ok_or(M::Error::custom("Cannot read t-component"))?;
        let x = access
            .next_element::<T>()?
            .ok_or(M::Error::custom("Cannot read x-component"))?;
        let y = access
            .next_element::<T>()?
            .ok_or(M::Error::custom("Cannot read y-component"))?;
        let z = access
            .next_element::<T>()?
            .ok_or(M::Error::custom("Cannot read z-component"))?;

        Ok(LorentzVector::from_args(t, x, y, z))
    }
}

impl<'de, T: Field + Deserialize<'de>> Deserialize<'de> for LorentzVector<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_any(LorentzVectorVisitor {
            _marker: PhantomData,
        })
    }
}
