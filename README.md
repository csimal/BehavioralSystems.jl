# BehavioralSystems

[![Build Status](https://github.com/csimal/BehavioralSystems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/csimal/BehavioralSystems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/csimal/BehavioralSystems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/csimal/BehavioralSystems.jl)

This package provides a collection of functions to analyze dynamical systems based on Behavioral Systems Theory. In particular, for LTI systems, we can use tools from linear algebra to determine characteristics of the system from a trajectory.

## Installation
This package is currently unregistered. To install it, enter the following command in the Julia REPL.
```
] add "https://github.com/csimal/BehavioralSystems.jl.git"
```

## Behavioral Approach to Dynamical Systems
Given a discrete-time LTI state space system $(A,B,C,D)$, with $n$ states, $m$ inputs and $p$ outputs, a *trajectory* of the system is a pair $w=(u,y) \in (\mathbb{R}^m)^\mathbb{N} \times (\mathbb{R}^p)^\mathbb{N}$, such that there is an $x\in (\mathbb{R}^n)^\mathbb{N}$ satisfying for all $t\in\mathbb{N}$,
- $x(t+1) = Ax(t) + Bu(t)$
- $y(t) = Cx(t) + Du(t)$

The set of all trajectories is called the *behavior* $\mathscr{B}$, and completely characterizes our dynamical system. Note that in many cases, the distinction between input and output is not so evident, and so trajectories are actually elements of $(\mathbb{R}^q)^\mathbb{N}$, with $q=m+p$.

For an LTI system, the behavior is a shift-invariant subspace of $(\mathbb{R}^q)^\mathbb{N}$, which allows us to study properties of the system from *restricted behaviors* $\mathscr{B}_T$ obtained from finite time-series of the system.

## Basic usage
See `/notebooks` for a showcase of the current features.