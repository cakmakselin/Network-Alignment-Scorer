# scoring_network_alignments
Programming assignment for "Fundamentals of Bioinformatics"

Task: In this assignment you are given three alignments of proteins of mus musculus ("mmu"), rattus norvegicus ("rno") and homo sapiens ("hsa"), which have been obtained through alignment of the corresponding protein-protein interaction networks. Since you have no clue about the performance of the used network alignment algorithm, you want to validate these alignments. One way of doing this is by computing a score based on shared GO terms: if the alignment is good then aligned protein pairs should have a lot of common GO terms.

Final script "score.py": Constructs mapping lists (each containing several dictionaries) by calling the function get_mapping on the two specified map files. Then we use these mapping lists as arguments to the function get_go_terms. The result of that function is the input of the compute_score function.

Examplary arguments via CLI: ./go_score.py alignments/rno-mmu.sif GO/rno.go GO/mmu.go mapping/rno.map mapping/mmu.map
