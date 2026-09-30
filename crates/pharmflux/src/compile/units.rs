//! Dimensional checking and conversion insertion. Runtime symbols remain in their declared units.
use super::*;
use pharmflux_core::units::{Dimension, Unit};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Number(Dimension),
    Boolean,
}
fn numeric(kind: Kind, path: &str) -> Result<Dimension, Error> {
    match kind {
        Kind::Number(d) => Ok(d),
        _ => Err(error(ErrorCode::Type, "expected numeric quantity", path)),
    }
}
fn same(a: Kind, b: Kind, path: &str) -> Result<(), Error> {
    if a == b {
        Ok(())
    } else {
        Err(error(ErrorCode::Unit, "expression dimensions differ", path))
    }
}
fn dimensionless(kind: Kind, path: &str) -> Result<(), Error> {
    same(kind, Kind::Number(Dimension::NONE), path)
}
fn scaled(expr: Expr, scale: f64) -> Expr {
    if scale == 1.0 {
        expr
    } else {
        Expr::binary(Binary::Multiply, Expr::number(scale), expr)
    }
}
fn lower(expr: &Expr, units: &BTreeMap<String, Unit>, path: &str) -> Result<(Expr, Kind), Error> {
    if let Expr::Literal { value } = expr {
        return Ok((
            expr.clone(),
            match value {
                Literal::Number(_) => Kind::Number(Dimension::NONE),
                Literal::Boolean(_) => Kind::Boolean,
            },
        ));
    }
    if let Expr::Symbol { name } = expr {
        let unit = units
            .get(name)
            .ok_or_else(|| error(ErrorCode::UnknownSymbol, "missing symbol unit", path))?;
        return Ok((
            scaled(expr.clone(), unit.scale.as_f64()),
            Kind::Number(unit.dimension),
        ));
    }
    let children = expr
        .children()
        .iter()
        .enumerate()
        .map(|(i, e)| lower(e, units, &format!("{path}.arguments[{i}]")))
        .collect::<Result<Vec<_>, _>>()?;
    let kinds: Vec<_> = children.iter().map(|(_, k)| *k).collect();
    let args: Vec<_> = children.into_iter().map(|(e, _)| e).collect();
    let (expr, kind) = match expr {
        Expr::Unary { operator, .. } => {
            let kind = if *operator == Unary::Not {
                same(kinds[0], Kind::Boolean, path)?;
                Kind::Boolean
            } else {
                Kind::Number(numeric(kinds[0], path)?)
            };
            (
                Expr::Unary {
                    operator: *operator,
                    arguments: Box::new([args[0].clone()]),
                },
                kind,
            )
        }
        Expr::Binary {
            operator,
            arguments,
        } => {
            let a = numeric(kinds[0], path)?;
            let b = numeric(kinds[1], path)?;
            let dimension = match operator {
                Binary::Add | Binary::Subtract | Binary::Modulo => {
                    same(kinds[0], kinds[1], path)?;
                    a
                }
                Binary::Multiply => a.multiply(b)?,
                Binary::Divide => a.divide(b)?,
                Binary::Power => {
                    dimensionless(kinds[1], path)?;
                    if a == Dimension::NONE {
                        a
                    } else {
                        let exponent = match &arguments[1] {
                            Expr::Literal {
                                value: Literal::Number(v),
                            } if v.is_finite()
                                && v.fract() == 0.0
                                && (-32.0..=32.0).contains(v) =>
                            {
                                *v as i16
                            }
                            _ => {
                                return Err(error(
                                    ErrorCode::Unit,
                                    "dimensional power requires an integer literal exponent",
                                    path,
                                ))
                            }
                        };
                        a.power(exponent)?
                    }
                }
            };
            (
                Expr::binary(*operator, args[0].clone(), args[1].clone()),
                Kind::Number(dimension),
            )
        }
        Expr::Compare { operator, .. } => {
            numeric(kinds[0], path)?;
            same(kinds[0], kinds[1], path)?;
            (
                Expr::Compare {
                    operator: *operator,
                    arguments: Box::new([args[0].clone(), args[1].clone()]),
                },
                Kind::Boolean,
            )
        }
        Expr::Boolean { operator, .. } => {
            for kind in &kinds {
                same(*kind, Kind::Boolean, path)?;
            }
            (
                Expr::Boolean {
                    operator: *operator,
                    arguments: args,
                },
                Kind::Boolean,
            )
        }
        Expr::Conditional { .. }
        | Expr::Call {
            name: Function::Ifelse,
            ..
        } => {
            same(kinds[0], Kind::Boolean, path)?;
            same(kinds[1], kinds[2], path)?;
            (
                Expr::conditional(args[0].clone(), args[1].clone(), args[2].clone()),
                kinds[1],
            )
        }
        Expr::Call { name, .. } => {
            let dimension = match name {
                Function::Abs => numeric(kinds[0], path)?,
                Function::Min | Function::Max => {
                    let dimension = numeric(kinds[0], path)?;
                    for kind in &kinds[1..] {
                        same(kinds[0], *kind, path)?;
                    }
                    dimension
                }
                Function::Sqrt => numeric(kinds[0], path)?.sqrt()?,
                Function::Hill => {
                    numeric(kinds[0], path)?;
                    same(kinds[0], kinds[2], path)?;
                    dimensionless(kinds[1], path)?;
                    Dimension::NONE
                }
                _ => {
                    dimensionless(kinds[0], path)?;
                    Dimension::NONE
                }
            };
            (Expr::call(*name, args), Kind::Number(dimension))
        }
        _ => unreachable!(),
    };
    Ok((expr, kind))
}
/// Validate expression shape/types and convert the resulting numeric value to the requested unit.
/// Inserted conversion work stays inside its original branch.
pub fn normalize(expr: &Expr, symbols: &[(Symbol, Unit)], output: Unit) -> Result<Expr, Error> {
    normalize_mode(expr, symbols, output, false)
}
pub(crate) fn normalize_expanded(
    expr: &Expr,
    symbols: &[(Symbol, Unit)],
    output: Unit,
) -> Result<Expr, Error> {
    normalize_mode(expr, symbols, output, true)
}
fn normalize_mode(
    expr: &Expr,
    symbols: &[(Symbol, Unit)],
    output: Unit,
    expanded: bool,
) -> Result<Expr, Error> {
    let names: Vec<_> = symbols.iter().map(|(symbol, _)| symbol.clone()).collect();
    if expanded {
        Program::compile_expanded(expr, &names)?;
    } else {
        Program::compile(expr, &names)?;
    }
    let units = symbols
        .iter()
        .map(|(symbol, unit)| (symbol.name.clone(), *unit))
        .collect();
    let (expr, kind) = lower(expr, &units, "root")?;
    same(kind, Kind::Number(output.dimension), "root")?;
    Ok(scaled(expr, 1.0 / output.scale.as_f64()))
}

/// Infer dimensions with unit scale one in the internal SI basis.
/// Explicit declared units may instead be selected by the caller.
pub fn infer(expr: &Expr, symbols: &[(Symbol, Unit)]) -> Result<Unit, Error> {
    infer_mode(expr, symbols, false)
}
pub(crate) fn infer_expanded(expr: &Expr, symbols: &[(Symbol, Unit)]) -> Result<Unit, Error> {
    infer_mode(expr, symbols, true)
}
fn infer_mode(expr: &Expr, symbols: &[(Symbol, Unit)], expanded: bool) -> Result<Unit, Error> {
    let names: Vec<_> = symbols.iter().map(|(s, _)| s.clone()).collect();
    if expanded {
        Program::compile_expanded(expr, &names)?;
    } else {
        Program::compile(expr, &names)?;
    }
    let units = symbols.iter().map(|(s, u)| (s.name.clone(), *u)).collect();
    let (_, kind) = lower(expr, &units, "root")?;
    Ok(Unit {
        dimension: numeric(kind, "root")?,
        scale: pharmflux_core::units::Rational::ONE,
    })
}
