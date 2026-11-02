use derive_more::Display;
use eyre::Report;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::node::GraphNodeKey;
use treetime_grid::{GridPointLimitExceeded, coordinates_coincide};
use treetime_utils::fmt::float::float_to_significant_digits;

pub(crate) trait WrapGridRequest<T> {
  fn wrap_grid_request(self, request: impl FnOnce() -> GridRequest) -> Result<T, Report>;
}

impl<T> WrapGridRequest<T> for Result<T, Report> {
  fn wrap_grid_request(self, request: impl FnOnce() -> GridRequest) -> Result<T, Report> {
    self.map_err(|report| {
      if report.downcast_ref::<GridRequest>().is_some() {
        report
      } else {
        report.wrap_err(request())
      }
    })
  }
}

#[derive(Clone, Debug, PartialEq, Display)]
#[display("When {}", operation.text(&key_label))]
pub(crate) struct GridRequest {
  operation: GridOperation,
  operands: Vec<GridOperand>,
}

impl GridRequest {
  pub(crate) fn new(operation: GridOperation, operands: impl IntoIterator<Item = GridOperand>) -> Self {
    Self {
      operation,
      operands: operands.into_iter().collect(),
    }
  }

  pub(crate) fn explain(
    &self,
    exceeded: &GridPointLimitExceeded,
    label: &dyn Fn(GraphNodeKey) -> String,
  ) -> GridExplanation {
    let (lo, hi) = exceeded.range;
    let dx = exceeded.dx;
    let grid = format!(
      "The grid spans [{}, {}] at spacing {}",
      number_text(lo),
      number_text(hi),
      number_text(dx)
    );
    if let GridOperation::BranchDistribution(_) = self.operation {
      return GridExplanation {
        operation: self.operation.text(label),
        attribution: format!("{grid}; the branch-length likelihood of the branch and the clock rate set both."),
        clock_rate_effect: Some(ClockRateEffect::Width),
      };
    }

    let spacing_setters = self
      .operands
      .iter()
      .filter(|operand| {
        operand
          .spacing
          .is_some_and(|spacing| coordinates_coincide(spacing, dx, spacing.max(dx)))
      })
      .collect_vec();
    let span_setters = self
      .operands
      .iter()
      .filter(|operand| {
        operand.range.is_some_and(|(start, end)| {
          coordinates_coincide(start, lo, start.abs().max(lo.abs()))
            && coordinates_coincide(end, hi, end.abs().max(hi.abs()))
        })
      })
      .collect_vec();

    let spacing = if spacing_setters.is_empty() {
      "no operand has this spacing".to_owned()
    } else {
      let setters = spacing_setters
        .iter()
        .map(|operand| operand.role.text_with_source(operand.origin.spacing, label))
        .join(" and ");
      format!("the spacing is that of {setters}")
    };
    let span = if span_setters.is_empty() {
      self.operation.kind().span_text().to_owned()
    } else {
      let setters = span_setters
        .iter()
        .map(|operand| operand.role.text_with_source(operand.origin.span, label))
        .join(" and ");
      format!("the span is that of {setters}")
    };
    let operands = self.operands.iter().map(|operand| operand.text(label)).join(", and ");

    let spacing_source = spacing_setters.iter().find_map(|operand| operand.origin.spacing);
    let span_source = span_setters
      .iter()
      .find_map(|operand| operand.origin.span)
      .or_else(|| self.operation.kind().combine(self.operands.iter().copied()).span);
    let clock_rate_effect = match (spacing_source, span_source) {
      (Some(GridSource::Branch(_)), Some(GridSource::Branch(_))) => Some(ClockRateEffect::Width),
      (Some(GridSource::Branch(_)), Some(GridSource::DateConstraint(_)) | None) => Some(ClockRateEffect::Narrow),
      (None, Some(GridSource::Branch(_))) => Some(ClockRateEffect::Wide),
      (Some(GridSource::DateConstraint(_)) | None, _) => None,
    };
    GridExplanation {
      operation: self.operation.text(label),
      attribution: format!("{grid}: {span}, and {spacing}. The operands are {operands}."),
      clock_rate_effect,
    }
  }
}

pub(crate) struct GridExplanation {
  pub operation: String,
  pub attribution: String,
  pub clock_rate_effect: Option<ClockRateEffect>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ClockRateEffect {
  Narrow,
  Wide,
  Width,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum GridOperation {
  BranchDistribution(GraphNodeKey),
  ChildMessages(GraphNodeKey),
  DateConstraint(GraphNodeKey),
  BackwardMessage(GraphNodeKey),
  ParentWithoutChild(GraphNodeKey),
  ForwardMessage(GraphNodeKey),
  Posterior(GraphNodeKey),
}

impl GridOperation {
  pub(crate) fn origin(self, operands: impl Iterator<Item = GridOperand> + Clone) -> GridOrigin {
    self.kind().combine(operands)
  }

  fn text(self, label: &dyn Fn(GraphNodeKey) -> String) -> String {
    match self {
      Self::BranchDistribution(node) => format!("building the time distribution of the branch above {}", label(node)),
      Self::ChildMessages(node) => format!("multiplying the time messages from the children of {}", label(node)),
      Self::DateConstraint(node) => format!("applying the date constraint of {}", label(node)),
      Self::BackwardMessage(node) => format!("sending the time message of {} backward to its parent", label(node)),
      Self::ParentWithoutChild(node) => format!(
        "removing the time message of {} from the time distribution of its parent",
        label(node)
      ),
      Self::ForwardMessage(node) => format!("sending the time message forward to {}", label(node)),
      Self::Posterior(node) => format!(
        "multiplying the forward time message of {} with its time distribution given its subtree",
        label(node)
      ),
    }
  }

  fn kind(self) -> OperationKind {
    match self {
      Self::BranchDistribution(_) => OperationKind::Construction,
      Self::BackwardMessage(_) | Self::ForwardMessage(_) => OperationKind::Convolution,
      Self::ChildMessages(_) | Self::DateConstraint(_) | Self::ParentWithoutChild(_) | Self::Posterior(_) => {
        OperationKind::Product
      },
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OperationKind {
  Construction,
  Convolution,
  Product,
}

impl OperationKind {
  fn combine(self, operands: impl Iterator<Item = GridOperand> + Clone) -> GridOrigin {
    let spacing = operands
      .clone()
      .filter_map(|operand| operand.spacing.map(|spacing| (spacing, operand.origin.spacing)))
      .min_by_key(|&(spacing, _)| OrderedFloat(spacing))
      .and_then(|(_, source)| source);
    let widths = operands.filter_map(|operand| operand.range.map(|(start, end)| (end - start, operand.origin.span)));
    let span = match self {
      Self::Convolution => widths.max_by_key(|&(width, _)| OrderedFloat(width)),
      Self::Product | Self::Construction => widths.min_by_key(|&(width, _)| OrderedFloat(width)),
    }
    .and_then(|(_, source)| source);
    GridOrigin { spacing, span }
  }

  fn span_text(self) -> &'static str {
    match self {
      Self::Convolution => "the span adds up the spans of the operands",
      Self::Product => "the span is the overlap of the operands",
      Self::Construction => "the span is that of the branch-length likelihood",
    }
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct GridOrigin {
  spacing: Option<GridSource>,
  span: Option<GridSource>,
}

impl GridOrigin {
  pub(crate) fn branch(node: GraphNodeKey) -> Self {
    Self::of(GridSource::Branch(node))
  }

  pub(crate) fn date_constraint(node: GraphNodeKey) -> Self {
    Self::of(GridSource::DateConstraint(node))
  }

  fn of(source: GridSource) -> Self {
    Self {
      spacing: Some(source),
      span: Some(source),
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GridSource {
  Branch(GraphNodeKey),
  DateConstraint(GraphNodeKey),
}

impl GridSource {
  fn role(self) -> OperandRole {
    match self {
      Self::Branch(node) => OperandRole::Branch(node),
      Self::DateConstraint(node) => OperandRole::DateConstraint(node),
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct GridOperand {
  role: OperandRole,
  range: Option<(f64, f64)>,
  spacing: Option<f64>,
  origin: GridOrigin,
}

impl GridOperand {
  pub(crate) fn of(role: OperandRole, distribution: &Distribution<NegLog>, origin: GridOrigin) -> Self {
    let spacing = match distribution {
      Distribution::Function(function) => Some(function.dx()),
      Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_) => None,
    };
    Self {
      role,
      range: distribution.time_bounds(),
      spacing,
      origin,
    }
  }

  fn text(&self, label: &dyn Fn(GraphNodeKey) -> String) -> String {
    let role = self.role.text(label);
    let range = self.range.map_or_else(String::new, |(start, end)| {
      format!(" over [{}, {}]", number_text(start), number_text(end))
    });
    let spacing = self
      .spacing
      .map_or_else(String::new, |spacing| format!(" at spacing {}", number_text(spacing)));
    format!("{role}{range}{spacing}")
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum OperandRole {
  Branch(GraphNodeKey),
  Subtree(GraphNodeKey),
  ChildMessage(GraphNodeKey),
  ChildProduct(GraphNodeKey),
  DateConstraint(GraphNodeKey),
  Posterior(GraphNodeKey),
  ParentWithoutChild(GraphNodeKey),
  ForwardMessage(GraphNodeKey),
}

impl OperandRole {
  fn text(self, label: &dyn Fn(GraphNodeKey) -> String) -> String {
    match self {
      Self::Branch(node) => format!("the time distribution of the branch above {}", label(node)),
      Self::Subtree(node) => format!("the time distribution of {} given its subtree", label(node)),
      Self::ChildMessage(node) => format!("the time message from {}", label(node)),
      Self::ChildProduct(node) => format!("the product of the time messages from the children of {}", label(node)),
      Self::DateConstraint(node) => format!("the date constraint of {}", label(node)),
      Self::Posterior(node) => format!("the time distribution of {}", label(node)),
      Self::ParentWithoutChild(node) => format!(
        "the time distribution of the parent of {} without the time message from {}",
        label(node),
        label(node)
      ),
      Self::ForwardMessage(node) => format!("the forward time message to {}", label(node)),
    }
  }

  fn text_with_source(self, source: Option<GridSource>, label: &dyn Fn(GraphNodeKey) -> String) -> String {
    let text = self.text(label);
    match source.map(GridSource::role) {
      Some(origin) if origin != self => format!("{text}, which takes it from {}", origin.text(label)),
      Some(_) | None => text,
    }
  }
}

fn key_label(key: GraphNodeKey) -> String {
  format!("node {key}")
}

fn number_text(value: f64) -> String {
  float_to_significant_digits(value, 6)
}
