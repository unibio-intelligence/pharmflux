//! Register bytecode. Branches are lowered as jumps; their domains stay guarded.
use pharmflux_core::expression::*;
use pharmflux_core::{Error, ErrorCode};
use std::collections::BTreeMap;
mod derivative;
pub mod units;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SymbolKind {
    Parameter,
    State,
    Time,
    Input,
}
#[derive(Clone, Debug)]
pub struct Symbol {
    pub name: String,
    pub kind: SymbolKind,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValueKind {
    Number,
    Boolean,
}
#[derive(Clone, Debug)]
enum Op {
    Constant(f64),
    Load(usize),
    Unary(Unary, usize),
    Binary(Binary, usize, usize),
    Compare(Compare, usize, usize),
    Call(Function, Vec<usize>),
    Copy(usize),
    JumpIfFalse(usize, usize),
    Jump(usize),
}
#[derive(Clone, Debug)]
struct Instruction {
    op: Op,
    destination: usize,
    path: String,
}
#[derive(Clone, Debug)]
pub struct Program {
    code: Vec<Instruction>,
    registers: usize,
    result: usize,
    inputs: usize,
    pub kind: ValueKind,
}
struct Compiler {
    code: Vec<Instruction>,
    next: usize,
    symbols: BTreeMap<String, (usize, SymbolKind)>,
    generated: bool,
}
#[derive(Clone, Copy)]
struct Register {
    index: usize,
    kind: ValueKind,
    dynamic: bool,
}
fn error(code: ErrorCode, message: impl Into<String>, path: &str) -> Error {
    let mut e = Error::new(code, message);
    e.expression = Some(path.into());
    e
}
fn check_size(expr: &Expr, depth: usize, nodes: &mut usize, limit: usize) -> Result<(), Error> {
    *nodes += 1;
    if depth > (if limit > 256 { 128 } else { 32 }) || *nodes > limit {
        return Err(error(
            ErrorCode::ExpressionLimit,
            "expression exceeds node or depth limit",
            "root",
        ));
    }
    for child in expr.children() {
        check_size(child, depth + 1, nodes, limit)?;
    }
    Ok(())
}
/// Output values are evaluated at a requested time/state, not integrated.
/// Preserve names, order and units while validating them as fixed readout inputs.
/// Never pass this view to dynamics, initial-state or dose compilation.
pub(crate) fn readout_symbols(
    symbols: &[(Symbol, pharmflux_core::units::Unit)],
) -> Vec<(Symbol, pharmflux_core::units::Unit)> {
    symbols
        .iter()
        .map(|(s, u)| {
            (
                Symbol {
                    name: s.name.clone(),
                    kind: SymbolKind::Parameter,
                },
                *u,
            )
        })
        .collect()
}
impl Program {
    pub fn compile(expr: &Expr, symbols: &[Symbol]) -> Result<Self, Error> {
        Self::compile_internal(expr, symbols, false)
    }
    pub(crate) fn compile_internal(
        expr: &Expr,
        symbols: &[Symbol],
        generated: bool,
    ) -> Result<Self, Error> {
        Self::compile_bounded(
            expr,
            symbols,
            generated,
            if generated { 16384 } else { 256 },
        )
    }
    /// Expanded user definitions retain dynamic-crossing restrictions while
    /// allowing bounded compiler expansion beyond one source expression.
    pub(crate) fn compile_expanded(expr: &Expr, symbols: &[Symbol]) -> Result<Self, Error> {
        Self::compile_bounded(expr, symbols, false, 16384)
    }
    fn compile_bounded(
        expr: &Expr,
        symbols: &[Symbol],
        generated: bool,
        limit: usize,
    ) -> Result<Self, Error> {
        check_size(expr, 1, &mut 0, limit)?;
        let mut c = Compiler {
            code: vec![],
            next: 0,
            symbols: BTreeMap::new(),
            generated,
        };
        for (index, symbol) in symbols.iter().enumerate() {
            if symbol.name.is_empty()
                || !symbol.name.bytes().enumerate().all(|(i, b)| {
                    b == b'_' || b.is_ascii_lowercase() || (i > 0 && b.is_ascii_digit())
                })
            {
                return Err(error(
                    ErrorCode::InvalidInput,
                    "expected ASCII snake_case symbol",
                    "root",
                ));
            }
            if c.symbols
                .insert(symbol.name.clone(), (index, symbol.kind))
                .is_some()
            {
                return Err(error(ErrorCode::InvalidInput, "duplicate symbol", "root"));
            }
        }
        let r = c.lower(expr, "root")?;
        Ok(Self {
            code: c.code,
            registers: c.next,
            result: r.index,
            inputs: symbols.len(),
            kind: r.kind,
        })
    }
    pub fn compile_with_units(
        expr: &Expr,
        symbols: &[(Symbol, pharmflux_core::units::Unit)],
        output: pharmflux_core::units::Unit,
    ) -> Result<Self, Error> {
        let normalized = units::normalize(expr, symbols, output)?;
        let names: Vec<_> = symbols.iter().map(|(symbol, _)| symbol.clone()).collect();
        // normalize validated the source. Its inserted conversions can exceed the
        // source node budget but cannot introduce a new switching condition.
        Self::compile_internal(&normalized, &names, true)
    }
    pub(crate) fn compile_expanded_with_units(
        expr: &Expr,
        symbols: &[(Symbol, pharmflux_core::units::Unit)],
        output: pharmflux_core::units::Unit,
    ) -> Result<Self, Error> {
        let normalized = units::normalize_expanded(expr, symbols, output)?;
        let names = symbols.iter().map(|(s, _)| s.clone()).collect::<Vec<_>>();
        Self::compile_internal(&normalized, &names, true)
    }
    pub fn derivative(
        expr: &Expr,
        symbols: &[Symbol],
        with_respect_to: &str,
    ) -> Result<Self, Error> {
        Self::derivative_internal(expr, symbols, with_respect_to, false)
    }
    pub(crate) fn derivative_internal(
        expr: &Expr,
        symbols: &[Symbol],
        with_respect_to: &str,
        normalized: bool,
    ) -> Result<Self, Error> {
        let original = Self::compile_internal(expr, symbols, normalized)?;
        if original.kind != ValueKind::Number {
            return Err(error(
                ErrorCode::Type,
                "cannot differentiate boolean output",
                "root",
            ));
        }
        if !symbols.iter().any(|s| s.name == with_respect_to) {
            return Err(error(
                ErrorCode::UnknownSymbol,
                "unknown differentiation symbol",
                "root",
            ));
        }
        // Retain primal domain checks even when derivative simplification yields zero.
        let checked = Expr::binary(
            Binary::Add,
            Expr::binary(Binary::Multiply, Expr::number(0.0), expr.clone()),
            derivative::differentiate(
                expr,
                with_respect_to,
                symbols
                    .iter()
                    .any(|s| s.name == with_respect_to && s.kind == SymbolKind::State),
            )?,
        );
        Self::compile_internal(&checked, symbols, true)
    }
    pub fn instruction_count(&self) -> usize {
        self.code.len()
    }
    pub(crate) fn workspace_size(&self) -> usize {
        self.registers
    }
    pub fn workspace(&self) -> Vec<f64> {
        vec![0.0; self.registers]
    }
    pub fn evaluate(&self, inputs: &[f64]) -> Result<f64, Error> {
        self.evaluate_into(inputs, &mut self.workspace())
    }
    pub fn evaluate_into(&self, inputs: &[f64], registers: &mut [f64]) -> Result<f64, Error> {
        if inputs.len() != self.inputs
            || registers.len() < self.registers
            || inputs.iter().any(|v| !v.is_finite())
        {
            return Err(error(
                ErrorCode::InvalidInput,
                "invalid input buffer or VM workspace",
                "root",
            ));
        }
        let mut pc = 0;
        while pc < self.code.len() {
            let instruction = &self.code[pc];
            let result = match &instruction.op {
                Op::Constant(v) => *v,
                Op::Load(i) => inputs[*i],
                Op::Copy(i) => registers[*i],
                Op::Jump(target) => {
                    pc = *target;
                    continue;
                }
                Op::JumpIfFalse(condition, target) => {
                    if registers[*condition] == 0.0 {
                        pc = *target;
                        continue;
                    }
                    pc += 1;
                    continue;
                }
                Op::Unary(op, i) => match op {
                    Unary::Positive => registers[*i],
                    Unary::Negative => -registers[*i],
                    Unary::Not => f64::from(registers[*i] == 0.0),
                },
                Op::Binary(op, a, b) => {
                    let (a, b) = (registers[*a], registers[*b]);
                    match op {
                        Binary::Add => a + b,
                        Binary::Subtract => a - b,
                        Binary::Multiply => a * b,
                        Binary::Divide => {
                            if b != 0.0 {
                                a / b
                            } else {
                                f64::NAN
                            }
                        }
                        Binary::Power => {
                            if (a >= 0.0 || b.fract() == 0.0) && !(a == 0.0 && b <= 0.0) {
                                libm::pow(a, b)
                            } else {
                                f64::NAN
                            }
                        }
                        Binary::Modulo => {
                            if b > 0.0 {
                                let r = a % b;
                                if r < 0.0 {
                                    let r = r + b;
                                    if r >= b {
                                        0.0
                                    } else {
                                        r
                                    }
                                } else {
                                    r
                                }
                            } else {
                                f64::NAN
                            }
                        }
                    }
                }
                Op::Compare(op, a, b) => {
                    let (a, b) = (registers[*a], registers[*b]);
                    f64::from(match op {
                        Compare::Lt => a < b,
                        Compare::Le => a <= b,
                        Compare::Gt => a > b,
                        Compare::Ge => a >= b,
                        Compare::Eq => a == b,
                        Compare::Ne => a != b,
                    })
                }
                Op::Call(f, args) => {
                    let x = registers[args[0]];
                    match f {
                        Function::Exp => libm::exp(x),
                        Function::Log => libm::log(x),
                        Function::Log10 => libm::log10(x),
                        Function::Sqrt => libm::sqrt(x),
                        Function::Sin => libm::sin(x),
                        Function::Cos => libm::cos(x),
                        Function::Tan => libm::tan(x),
                        Function::Tanh => libm::tanh(x),
                        Function::Abs => x.abs(),
                        Function::Floor => libm::floor(x),
                        Function::Ceil => libm::ceil(x),
                        Function::Min => args.iter().fold(x, |v, i| v.min(registers[*i])),
                        Function::Max => args.iter().fold(x, |v, i| v.max(registers[*i])),
                        Function::Hill => {
                            let n = registers[args[1]];
                            let k = registers[args[2]];
                            if n <= 0.0 || k <= 0.0 {
                                f64::NAN
                            } else if x <= 0.0 {
                                0.0
                            } else {
                                let z = n * (libm::log(x) - libm::log(k));
                                if z >= 0.0 {
                                    1.0 / (1.0 + libm::exp(-z))
                                } else {
                                    let q = libm::exp(z);
                                    q / (1.0 + q)
                                }
                            }
                        }
                        Function::Ifelse => unreachable!("ifelse is lowered into branches"),
                    }
                }
            };
            if !result.is_finite() {
                return Err(error(
                    ErrorCode::Domain,
                    format!("non-finite or out-of-domain {:?}", instruction.op),
                    &instruction.path,
                ));
            }
            registers[instruction.destination] = result;
            pc += 1;
        }
        Ok(registers[self.result])
    }
}
impl Compiler {
    fn emit(&mut self, op: Op, kind: ValueKind, dynamic: bool, path: &str) -> Register {
        let index = self.next;
        self.next += 1;
        self.code.push(Instruction {
            op,
            destination: index,
            path: path.into(),
        });
        Register {
            index,
            kind,
            dynamic,
        }
    }
    fn require(&self, r: Register, kind: ValueKind, path: &str) -> Result<(), Error> {
        if r.kind != kind {
            Err(error(
                ErrorCode::Type,
                "operator operand type mismatch",
                path,
            ))
        } else {
            Ok(())
        }
    }
    fn branch(
        &mut self,
        condition: Register,
        yes: &Expr,
        no: &Expr,
        path: &str,
    ) -> Result<Register, Error> {
        self.require(condition, ValueKind::Boolean, path)?;
        if condition.dynamic && !self.generated {
            return Err(error(
                ErrorCode::Unsupported,
                "state/time/input-dependent switches require a qualified event contract",
                path,
            ));
        }
        let test = self.code.len();
        self.emit(
            Op::JumpIfFalse(condition.index, 0),
            ValueKind::Boolean,
            false,
            path,
        );
        let y = self.lower(yes, &format!("{path}.then"))?;
        let result = self.emit(Op::Copy(y.index), y.kind, y.dynamic, path);
        let jump = self.code.len();
        self.emit(Op::Jump(0), ValueKind::Boolean, false, path);
        let alternative = self.code.len();
        self.code[test].op = Op::JumpIfFalse(condition.index, alternative);
        let n = self.lower(no, &format!("{path}.else"))?;
        self.require(n, y.kind, path)?;
        self.code.push(Instruction {
            op: Op::Copy(n.index),
            destination: result.index,
            path: path.into(),
        });
        self.code[jump].op = Op::Jump(self.code.len());
        Ok(Register {
            dynamic: condition.dynamic || y.dynamic || n.dynamic,
            ..result
        })
    }
    fn lower(&mut self, expr: &Expr, path: &str) -> Result<Register, Error> {
        use ValueKind::{Boolean as Bool, Number};
        match expr {
            Expr::Literal { value } => match value {
                Literal::Number(v) => {
                    if v.is_finite() {
                        Ok(self.emit(Op::Constant(*v), Number, false, path))
                    } else {
                        Err(error(ErrorCode::Domain, "non-finite literal", path))
                    }
                }
                Literal::Boolean(v) => {
                    Ok(self.emit(Op::Constant(f64::from(*v)), Bool, false, path))
                }
            },
            Expr::Symbol { name } => {
                let (index, kind) = *self.symbols.get(name).ok_or_else(|| {
                    error(
                        ErrorCode::UnknownSymbol,
                        format!("unknown symbol {name}"),
                        path,
                    )
                })?;
                Ok(self.emit(Op::Load(index), Number, kind != SymbolKind::Parameter, path))
            }
            Expr::Conditional { arguments } => {
                let c = self.lower(&arguments[0], &format!("{path}.condition"))?;
                self.branch(c, &arguments[1], &arguments[2], path)
            }
            Expr::Call {
                name: Function::Ifelse,
                arguments,
            } => {
                if arguments.len() != 3 {
                    return Err(error(
                        ErrorCode::InvalidInput,
                        "ifelse requires three arguments",
                        path,
                    ));
                }
                let c = self.lower(&arguments[0], &format!("{path}.condition"))?;
                self.branch(c, &arguments[1], &arguments[2], path)
            }
            Expr::Boolean {
                operator,
                arguments,
            } => {
                if arguments.len() < 2 {
                    return Err(error(
                        ErrorCode::InvalidInput,
                        "boolean operator requires at least two arguments",
                        path,
                    ));
                }
                let mut expression = arguments.last().unwrap().clone();
                for arg in arguments[..arguments.len() - 1].iter().rev() {
                    expression = match operator {
                        Boolean::And => Expr::conditional(
                            arg.clone(),
                            expression,
                            Expr::Literal {
                                value: Literal::Boolean(false),
                            },
                        ),
                        Boolean::Or => Expr::conditional(
                            arg.clone(),
                            Expr::Literal {
                                value: Literal::Boolean(true),
                            },
                            expression,
                        ),
                    };
                }
                let r = self.lower(&expression, path)?;
                self.require(r, Bool, path)?;
                Ok(r)
            }
            _ => {
                let args: Vec<_> = expr
                    .children()
                    .iter()
                    .enumerate()
                    .map(|(i, e)| self.lower(e, &format!("{path}.arguments[{i}]")))
                    .collect::<Result<_, _>>()?;
                let dynamic = args.iter().any(|r| r.dynamic);
                let (op, kind) = match expr {
                    Expr::Unary { operator, .. } => {
                        self.require(
                            args[0],
                            if *operator == Unary::Not {
                                Bool
                            } else {
                                Number
                            },
                            path,
                        )?;
                        (
                            Op::Unary(*operator, args[0].index),
                            if *operator == Unary::Not {
                                Bool
                            } else {
                                Number
                            },
                        )
                    }
                    Expr::Binary { operator, .. } => {
                        for r in &args {
                            self.require(*r, Number, path)?;
                        }
                        if *operator == Binary::Modulo && dynamic && !self.generated {
                            return Err(error(
                                ErrorCode::Unsupported,
                                "dynamic modulo requires declared boundaries",
                                path,
                            ));
                        }
                        (Op::Binary(*operator, args[0].index, args[1].index), Number)
                    }
                    Expr::Compare { operator, .. } => {
                        for r in &args {
                            self.require(*r, Number, path)?;
                        }
                        (Op::Compare(*operator, args[0].index, args[1].index), Bool)
                    }
                    Expr::Call { name, .. } => {
                        let valid = match name {
                            Function::Hill => args.len() == 3,
                            Function::Min | Function::Max => args.len() >= 2,
                            _ => args.len() == 1,
                        };
                        if !valid {
                            return Err(error(
                                ErrorCode::InvalidInput,
                                "invalid function arity",
                                path,
                            ));
                        }
                        for r in &args {
                            self.require(*r, Number, path)?;
                        }
                        if matches!(name, Function::Floor | Function::Ceil)
                            && dynamic
                            && !self.generated
                        {
                            return Err(error(
                                ErrorCode::Unsupported,
                                "dynamic rounding requires declared boundaries",
                                path,
                            ));
                        }
                        (
                            Op::Call(*name, args.iter().map(|r| r.index).collect()),
                            Number,
                        )
                    }
                    _ => unreachable!(),
                };
                Ok(self.emit(op, kind, dynamic, path))
            }
        }
    }
}

/// Internal expression form for augmented sensitivity equations. The caller
/// validates the original model first; retain its domain checks in each partial.
pub(crate) fn derivative_expression(expr: &Expr, name: &str) -> Result<Expr, Error> {
    Ok(Expr::binary(
        Binary::Add,
        Expr::binary(Binary::Multiply, Expr::number(0.), expr.clone()),
        derivative::differentiate(expr, name, false)?,
    ))
}
