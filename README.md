# Project 5
Needleman Wunsch Algorithm


# Assignment Overview
The purpose of this assigment is to have you implement the Needleman-Wunsch global pairwise sequence alignment algorithm (dynamic programming).
See this [video](https://www.youtube.com/watch?v=NqYY0PJbD3s) for a walk through of the algorithm implementation. In addition, this is also a helpful resource: https://www.bioinformaticsalgorithms.org/bioinformatics-chapter-5.

# Assignment Tasks
## Coding Assessment
**Note: All modules you need have already been imported.**
* [X] Complete the `NeedlemanWunsch.align` method found in the align/align.py 
	* Finish the method for filling in the alignment and gap matrices for Needleman-Wunsch.
* [X] Complete the `NeedlemanWunsch._backtrace` method found in align/align.py
  * Use the heuristic you have chosen to implement in the the `NeedlemanWunsch.align` method to implement the backtracing procedure.
* [X] Complete the `main` function in main.py to 
    1. align all provided species BRD2 sequences to the human BRD2 sequence and print the species in order of most similar to least similar with respect to human BRD2.
    2. print the alignment scores corresponding to each species alignment to the human BRD2 sequence.

## Software Development Assessment
### Unit Tests
* [X] Complete the `test_nw_alignment` function in test/test_align.py to test for proper matrix filling in your `NeedlemanWunsch.align` method.
* [X] Complete the `test_nw_backtrace` function in test/test_align.py to test for proper backtracing in your `NeedlemanWunsch._backtrace` method.

Note: To check that you have correclty implemented your algorithm the following information is provided
* test_seq3.fa and test_seq4.fa should have an alignment score of **17** and an alignment of:

	MAVHQLIRRP
	
	M---QLIRHP


[X] Ensure that it runs pytest

### Pip Installable
* [X] make .toml file with flit and ensure that your package can be installed with pip

# Getting Started
To get started you will need to fork this repository onto your own Github account. Work on the codebase from your own repo and commit changes. 

The following packages will be needed:
* numpy
* pytest

# Completing the assignment
Make sure to push all your code to Github, ensure that your unit tests are correct, and submit a link to your Github through the Google classroom assignment.

# Grading
## Code (6 points)
* Pairwise global alignment works properly (6)
    * Correct implementation of Needleman-Wunsch algorithm (4)
    * Produces correct order of species in main.py (1) 
    * Produces correct NW alignment scores in main.py (1)

## Unit tests (3 points)
* `test_nw_alignment` function properly checks that matrices are filled in correctly for alignment of test_seq1.fa and test_seq2.fa (1)
* `test_nw_backtrace` function properly checks that backtrace works correctly (1)
* Ensure functionality with pytest (1)
## Style (1 points)
* Readable code with clear comments and method descriptions (1)
## Extra credit (0.5)
* Github actions/workflow (0.5)
