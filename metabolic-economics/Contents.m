% MATLAB code for Metabolic Economics
%
% The toolbox is based on the Metabolic Network Toolbox and uses
% its functions and data structures. Main purposes are:
%
%   o Compute economical flux distributions
%   o Test if flux distributions are economical (and thermo-feasible)
%   o Compute economic potentials in agreement with a given flux distribution
%   o Realise an economical flux distribution by an enzyme-balanced kinetic model
%   o Optimise the enzyme levels in a kinetic model, leading to an enzyme-optimal state
%
% Options for all functions are stored in common
% matlab structs called 'cba_options' and 'cba_constraints'.
%
%
% Main functions (in directory 'cba')
%
% Options
%   cba_default_options     - Set options in data structures 'cba_options' and 'cba_constraints'
%
% Feasible flux distribution 
%   cba_adjust_fluxes       - Correct a flux distribution by removing economically infeasible cycles
%
% Check flux distribution (based on different types of information given)
%   cba_feasible_cycle      - Check economical flux distribution (non-beneficial cycle criterion)
%   cba_feasible_efm        - Check economical flux distribution (given non-beneficial flux modes)
%   cba_feasible_lp         - Check flux mode for EFA criteria; choose enzyme costs and thermodynamic forces
%
% Compute economic potentials given fluxes
%   cba_homogeneous_cost    - Determine economic potentials from the principle of even enzyme cost
%
% Compute fluxes and economic potentials
%   cba_economic_state            - Find fluxes, chemical and economic potentials satisfying FBA/EBA/CBA constraints
%   cba_check               - Check flux mode, chem. pot. differences, and spec. flux costs for feasibility
%
% Enzyme-balanced kinetic models
%   cba_reconstruct_model   - Reconstruct enzyme-balanced model from economical flux mode
%   cba_economic_potentials - Compute economic potentials and their sensitivities in a kinetic model
%   cba_check_model         - Check kinetic model for being in an enzyme-optimal state
%
% Further utility functions can be found in 'cba_utils'
%
%
% MATLAB Toolbox required:
%   Metabolic Network Toolbox (https://github.com/wolframliebermeister/metabolic-network-toolbox)
%   SBMLtoolbox    - SBML import / export  (see http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox  - SBtab format (https://github.com/wolframliebermeister/sbtab-matlab)
%
% Copyright (C) 2014
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>
