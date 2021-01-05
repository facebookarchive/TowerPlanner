# General
Tower Planner is a light-weight, structural analysis and optimization toolkit designed for cell towers. It provides rapid structural capacity evaluation, tower deformation due to wind-loading and generates cost-optimal tower designs given high-level requirements. 

A web-app version is also available on www.fbctower.info

## Project Overview
Tower structural analysis is an integral part of evaluating a tower asset for multi-tenancy potential and site-valuation. This activity is often time-consuming and expensive requiring accurate design details, drawings and significant structural engineering expertise. Tower Planner provides rapid analysis (in a few mins) that could help assist tower companies make faster business decisions towards procurement and leasing with miminal, basic data on the tower asset. Additionally, Tower Planner also helps design low-cost towers for asset-price sensitive scenarios (small cells and supercells) given basic, high-level requirements. 

## References
1. Aeroelastic Preliminary-Design Optimization of Communication Tower Structures, https://research.fb.com/publications/aeroelastic-preliminary-design-optimization-of-communication-tower-structures/
2. Frame3DD, http://frame3dd.sourceforge.net/
3. Midaco, http://www.midaco-solver.com/
4. Supercell, https://engineering.fb.com/2020/12/03/connectivity/supercell-reaching-new-heights-for-wider-connectivity/

## Requirements
Tower Planner was developed using MATLAB 2020a and requires MATLAB v2020a or later. Tower Planner also depends on other libraries (Frame3DD, Midaco) that are provided with this distribution. 

## Description
Tower Planner consists of routines to build structural models in Frame3DD [2], an open-source truss solver, from basic tower geometry data such as tower height, type of tower (self-supported, monopole or guyed), etc. Nodes in 3D space are first located, using which truss elements are constructed. Structural properties (stiffness and strength) and wind loading are assigned to each element. Structural properties and loading are specified in close alignment to the structural standard for cell towers: TIA-222-G. The finite-element model is solved for from which displacements and stresses developed in the structure are calculated. This information is useful to understand structural capacity available in the tower for additional loading. Furthermore, design constraints such as stress margins, guy wire slack or angular displacements of antennas exceeding allowable limits may be evaluated for the purpose of generating optimal designs. Optimization is performed with Midaco [3] using a discrete, heurestics-based method. Seperate optimization-cases are specified for monopoles, self-supported and guyed towers. More information is available in Ref. [1]. 

## Examples
Try Tower Planner from tests available in the /tests folder. Default values are specified using the app = loadDefaults() command. Values may then be modified as desired. Run results = execute(app, mode) to perform an analysis or design a tower depending on the mode variable ('analysis' or 'design'). 


## Join the TowerPlanner community
See the CONTRIBUTING file for how to help out.

## License
By contributing to TowerPlanner, you agree that your contributions will be licensed
under the LICENSE or COPYING file in the root directory of this source tree.
