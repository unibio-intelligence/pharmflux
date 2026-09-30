//! Explicit case-sensitive UCUM subset. Scales are exact positive rationals.
use crate::{Error, ErrorCode};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Rational {
    numerator: u128,
    denominator: u128,
}
fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}
fn failure(code: ErrorCode, message: &str) -> Error {
    Error::new(code, message)
}
impl Rational {
    pub fn new(numerator: u128, denominator: u128) -> Result<Self, Error> {
        if numerator == 0 || denominator == 0 {
            return Err(failure(ErrorCode::Unit, "unit scale must be positive"));
        }
        let divisor = gcd(numerator, denominator);
        Ok(Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        })
    }
    pub const ONE: Self = Self {
        numerator: 1,
        denominator: 1,
    };
    pub fn numerator(self) -> u128 {
        self.numerator
    }
    pub fn denominator(self) -> u128 {
        self.denominator
    }
    pub fn as_f64(self) -> f64 {
        self.numerator as f64 / self.denominator as f64
    }
    pub fn multiply(self, other: Self) -> Result<Self, Error> {
        let a = gcd(self.numerator, other.denominator);
        let b = gcd(other.numerator, self.denominator);
        let numerator = (self.numerator / a).checked_mul(other.numerator / b);
        let denominator = (self.denominator / b).checked_mul(other.denominator / a);
        Self::new(
            numerator.ok_or_else(|| failure(ErrorCode::Unit, "unit scale overflow"))?,
            denominator.ok_or_else(|| failure(ErrorCode::Unit, "unit scale overflow"))?,
        )
    }
    pub fn divide(self, other: Self) -> Result<Self, Error> {
        self.multiply(Self {
            numerator: other.denominator,
            denominator: other.numerator,
        })
    }
    pub fn power(self, exponent: i16) -> Result<Self, Error> {
        if !(-32..=32).contains(&exponent) {
            return Err(failure(ErrorCode::Unit, "unit exponent exceeds 32"));
        }
        let mut out = Self::ONE;
        let base = if exponent < 0 {
            Self {
                numerator: self.denominator,
                denominator: self.numerator,
            }
        } else {
            self
        };
        for _ in 0..exponent.unsigned_abs() {
            out = out.multiply(base)?;
        }
        Ok(out)
    }
}
/// SI kg, m, s, mol, K, A, cd exponents plus an opaque international-activity
/// exponent. The last coordinate is an internal compatibility guard, not an SI
/// definition of IU. Matching activity quantities must refer to the same assay
/// and substance; there is no automatic activity-to-mass/amount conversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Dimension(pub [i16; 8]);
impl Dimension {
    pub const NONE: Self = Self([0; 8]);
    pub fn multiply(self, other: Self) -> Result<Self, Error> {
        let mut out = [0; 8];
        for (i, value) in out.iter_mut().enumerate() {
            *value = self.0[i]
                .checked_add(other.0[i])
                .filter(|n| (-128..=128).contains(n))
                .ok_or_else(|| failure(ErrorCode::Unit, "dimension exponent overflow"))?;
        }
        Ok(Self(out))
    }
    pub fn power(self, exponent: i16) -> Result<Self, Error> {
        let mut out = [0; 8];
        for (i, value) in out.iter_mut().enumerate() {
            *value = self.0[i]
                .checked_mul(exponent)
                .filter(|n| (-128..=128).contains(n))
                .ok_or_else(|| failure(ErrorCode::Unit, "dimension exponent overflow"))?;
        }
        Ok(Self(out))
    }
    pub fn divide(self, other: Self) -> Result<Self, Error> {
        self.multiply(other.power(-1)?)
    }
    pub fn sqrt(self) -> Result<Self, Error> {
        if self.0.iter().any(|n| n % 2 != 0) {
            Err(failure(
                ErrorCode::Unit,
                "sqrt requires even dimension exponents",
            ))
        } else {
            Ok(Self(self.0.map(|n| n / 2)))
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Unit {
    pub dimension: Dimension,
    pub scale: Rational,
}
impl Unit {
    pub const ONE: Self = Self {
        dimension: Dimension::NONE,
        scale: Rational::ONE,
    };
    pub fn parse(code: &str) -> Result<Self, Error> {
        if code.is_empty() || code.len() > 256 || !code.is_ascii() {
            return Err(failure(
                ErrorCode::UnknownUnit,
                "invalid or oversized unit code",
            ));
        }
        let mut parser = Parser {
            bytes: code.as_bytes(),
            position: 0,
        };
        let unit = parser.term(0)?;
        if parser.position != parser.bytes.len() {
            return Err(failure(ErrorCode::UnknownUnit, "unrecognized unit suffix"));
        }
        Ok(unit)
    }
    pub fn multiply(self, other: Self) -> Result<Self, Error> {
        Ok(Self {
            dimension: self.dimension.multiply(other.dimension)?,
            scale: self.scale.multiply(other.scale)?,
        })
    }
    pub fn divide(self, other: Self) -> Result<Self, Error> {
        Ok(Self {
            dimension: self.dimension.divide(other.dimension)?,
            scale: self.scale.divide(other.scale)?,
        })
    }
    pub fn power(self, exponent: i16) -> Result<Self, Error> {
        Ok(Self {
            dimension: self.dimension.power(exponent)?,
            scale: self.scale.power(exponent)?,
        })
    }
    pub fn conversion_to(self, target: Self) -> Result<Rational, Error> {
        if self.dimension != target.dimension {
            return Err(failure(
                ErrorCode::Unit,
                "incompatible dimensions; mass/amount conversion requires explicit molar mass",
            ));
        }
        self.scale.divide(target.scale)
    }
    pub fn convert(self, value: f64, target: Self) -> Result<f64, Error> {
        if !value.is_finite() {
            return Err(failure(ErrorCode::InvalidInput, "non-finite quantity"));
        }
        let result = value * self.conversion_to(target)?.as_f64();
        if !result.is_finite() {
            Err(failure(ErrorCode::Domain, "quantity conversion overflow"))
        } else {
            Ok(result)
        }
    }
}
fn atom(code: &str) -> Result<Unit, Error> {
    if code == "[iU]" || code == "[IU]" {
        return Ok(Unit {
            dimension: Dimension([0, 0, 0, 0, 0, 0, 0, 1]),
            scale: Rational::ONE,
        });
    }
    let (dimension, numerator, denominator) = match code {
        "1" => ([0; 7], 1, 1),
        "%" => ([0; 7], 1, 100),
        "kg" => ([1, 0, 0, 0, 0, 0, 0], 1, 1),
        "g" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000),
        "mg" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000000),
        "ug" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000000000),
        "ng" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000000000000),
        "pg" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000000000000000),
        "fg" => ([1, 0, 0, 0, 0, 0, 0], 1, 1000000000000000000),
        "m" => ([0, 1, 0, 0, 0, 0, 0], 1, 1),
        "cm" => ([0, 1, 0, 0, 0, 0, 0], 1, 100),
        "mm" => ([0, 1, 0, 0, 0, 0, 0], 1, 1000),
        "um" => ([0, 1, 0, 0, 0, 0, 0], 1, 1000000),
        "nm" => ([0, 1, 0, 0, 0, 0, 0], 1, 1000000000),
        "L" | "l" => ([0, 3, 0, 0, 0, 0, 0], 1, 1000),
        "dL" => ([0, 3, 0, 0, 0, 0, 0], 1, 10000),
        "mL" => ([0, 3, 0, 0, 0, 0, 0], 1, 1000000),
        "uL" => ([0, 3, 0, 0, 0, 0, 0], 1, 1000000000),
        "nL" => ([0, 3, 0, 0, 0, 0, 0], 1, 1000000000000),
        "pL" => ([0, 3, 0, 0, 0, 0, 0], 1, 1000000000000000),
        "s" => ([0, 0, 1, 0, 0, 0, 0], 1, 1),
        "ms" => ([0, 0, 1, 0, 0, 0, 0], 1, 1000),
        "us" => ([0, 0, 1, 0, 0, 0, 0], 1, 1000000),
        "ns" => ([0, 0, 1, 0, 0, 0, 0], 1, 1000000000),
        "min" => ([0, 0, 1, 0, 0, 0, 0], 60, 1),
        "h" => ([0, 0, 1, 0, 0, 0, 0], 3600, 1),
        "d" => ([0, 0, 1, 0, 0, 0, 0], 86400, 1),
        "wk" => ([0, 0, 1, 0, 0, 0, 0], 604800, 1),
        "a" => ([0, 0, 1, 0, 0, 0, 0], 31557600, 1),
        "mol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1),
        "mmol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1000),
        "umol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1000000),
        "nmol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1000000000),
        "pmol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1000000000000),
        "fmol" => ([0, 0, 0, 1, 0, 0, 0], 1, 1000000000000000),
        "K" => ([0, 0, 0, 0, 1, 0, 0], 1, 1),
        "A" => ([0, 0, 0, 0, 0, 1, 0], 1, 1),
        "cd" => ([0, 0, 0, 0, 0, 0, 1], 1, 1),
        "rad" | "sr" => ([0; 7], 1, 1),
        "Hz" | "Bq" => ([0, 0, -1, 0, 0, 0, 0], 1, 1),
        "N" => ([1, 1, -2, 0, 0, 0, 0], 1, 1),
        "Pa" => ([1, -1, -2, 0, 0, 0, 0], 1, 1),
        "J" => ([1, 2, -2, 0, 0, 0, 0], 1, 1),
        "W" => ([1, 2, -3, 0, 0, 0, 0], 1, 1),
        "C" => ([0, 0, 1, 0, 0, 1, 0], 1, 1),
        "V" => ([1, 2, -3, 0, 0, -1, 0], 1, 1),
        "Ohm" => ([1, 2, -3, 0, 0, -2, 0], 1, 1),
        "S" => ([-1, -2, 3, 0, 0, 2, 0], 1, 1),
        "F" => ([-1, -2, 4, 0, 0, 2, 0], 1, 1),
        "Wb" => ([1, 2, -2, 0, 0, -1, 0], 1, 1),
        "T" => ([1, 0, -2, 0, 0, -1, 0], 1, 1),
        "H" => ([1, 2, -2, 0, 0, -2, 0], 1, 1),
        "lm" => ([0, 0, 0, 0, 0, 0, 1], 1, 1),
        "lx" => ([0, -2, 0, 0, 0, 0, 1], 1, 1),
        "Gy" | "Sv" => ([0, 2, -2, 0, 0, 0, 0], 1, 1),
        "kat" => ([0, 0, -1, 1, 0, 0, 0], 1, 1),
        _ => {
            return Err(failure(
                ErrorCode::UnknownUnit,
                "unit atom is outside the supported UCUM subset",
            ))
        }
    };
    Ok(Unit {
        dimension: Dimension([
            dimension[0],
            dimension[1],
            dimension[2],
            dimension[3],
            dimension[4],
            dimension[5],
            dimension[6],
            0,
        ]),
        scale: Rational::new(numerator, denominator)?,
    })
}
struct Parser<'a> {
    bytes: &'a [u8],
    position: usize,
}
impl Parser<'_> {
    fn term(&mut self, depth: usize) -> Result<Unit, Error> {
        if depth > 16 {
            return Err(failure(ErrorCode::UnknownUnit, "unit nesting exceeds 16"));
        }
        let mut value = self.factor(depth)?;
        while let Some(operator @ (b'.' | b'/')) = self.bytes.get(self.position).copied() {
            self.position += 1;
            let rhs = self.factor(depth)?;
            value = if operator == b'.' {
                value.multiply(rhs)?
            } else {
                value.divide(rhs)?
            };
        }
        Ok(value)
    }
    fn factor(&mut self, depth: usize) -> Result<Unit, Error> {
        let mut value = if self.bytes.get(self.position) == Some(&b'(') {
            self.position += 1;
            let value = self.term(depth + 1)?;
            if self.bytes.get(self.position) != Some(&b')') {
                return Err(failure(ErrorCode::UnknownUnit, "unclosed unit group"));
            }
            self.position += 1;
            value
        } else if self.bytes.get(self.position) == Some(&b'[') {
            let start = self.position;
            while self.bytes.get(self.position).is_some_and(|c| *c != b']') {
                self.position += 1;
            }
            if self.bytes.get(self.position) != Some(&b']') {
                return Err(failure(ErrorCode::UnknownUnit, "unclosed bracketed unit"));
            }
            self.position += 1;
            atom(
                std::str::from_utf8(&self.bytes[start..self.position])
                    .map_err(|_| failure(ErrorCode::UnknownUnit, "non-ASCII unit"))?,
            )?
        } else {
            let start = self.position;
            while self
                .bytes
                .get(self.position)
                .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'%')
            {
                self.position += 1;
            }
            if start == self.position
                && self
                    .bytes
                    .get(self.position)
                    .is_some_and(u8::is_ascii_digit)
            {
                while self
                    .bytes
                    .get(self.position)
                    .is_some_and(u8::is_ascii_digit)
                {
                    self.position += 1;
                }
                let numerator = std::str::from_utf8(&self.bytes[start..self.position])
                    .ok()
                    .and_then(|s| s.parse::<u128>().ok())
                    .ok_or_else(|| failure(ErrorCode::Unit, "unit scale overflow"))?;
                // Consume the entire integer as a factor, never as an exponent
                // of 1. Bare numeric powers and implicit products stay excluded;
                // use a group for a power, e.g. (1000)-1.
                return Ok(Unit {
                    dimension: Dimension::NONE,
                    scale: Rational::new(numerator, 1)?,
                });
            }
            if start == self.position {
                return Err(failure(ErrorCode::UnknownUnit, "expected unit atom"));
            }
            atom(
                std::str::from_utf8(&self.bytes[start..self.position])
                    .map_err(|_| failure(ErrorCode::UnknownUnit, "non-ASCII unit"))?,
            )?
        };
        let start = self.position;
        if self
            .bytes
            .get(self.position)
            .is_some_and(|c| *c == b'-' || *c == b'+')
        {
            self.position += 1;
        }
        let digits = self.position;
        while self
            .bytes
            .get(self.position)
            .is_some_and(u8::is_ascii_digit)
        {
            self.position += 1;
        }
        if self.position > digits {
            let exponent = std::str::from_utf8(&self.bytes[start..self.position])
                .ok()
                .and_then(|s| s.parse::<i16>().ok())
                .ok_or_else(|| failure(ErrorCode::Unit, "invalid unit exponent"))?;
            value = value.power(exponent)?;
        } else if self.position != start {
            return Err(failure(
                ErrorCode::UnknownUnit,
                "unit exponent has no digits",
            ));
        }
        Ok(value)
    }
}
