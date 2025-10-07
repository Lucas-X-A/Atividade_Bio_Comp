import numpy as np
from needleman_wunsch import needleman_wunsch

def traceback_from_matrix(seq1, seq2, score_matrix, match_score=1, mismatch_score=-1, gap_penalty=-2):
    """
    Obtém o alinhamento global a partir de uma matriz de programação dinâmica pré-calculada.

    Args:
        seq1 (str): A primeira sequência original.
        seq2 (str): A segunda sequência original.
        score_matrix (numpy.ndarray): A matriz de pontuação gerada pelo Needleman-Wunsch.
        match_score (int): Pontuação para um match.
        mismatch_score (int): Pontuação para um mismatch.
        gap_penalty (int): Pontuação para um gap (lacuna).

    Returns:
        tuple: Uma tupla contendo a primeira sequência alinhada (str) e
               a segunda sequência alinhada (str).
    """
    # Inicializa as sequências alinhadas como strings vazias
    
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    # Começa do canto inferior direito da matriz
    i, j = len(seq2), len(seq1)

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        # Verifica se o caminho veio da diagonal (match/mismatch)
        # Este movimento é preferencial em caso de empate de pontuação.
        if i > 0 and j > 0 and current_score == score_matrix[i-1][j-1] + (match_score if seq1[j-1] == seq2[i-1] else mismatch_score):
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            i -= 1
            j -= 1
        # Verifica se o caminho veio de cima (gap em seq1)
        elif i > 0 and current_score == score_matrix[i-1][j] + gap_penalty:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            i -= 1
        # Verifica se o caminho veio da esquerda (gap em seq2)
        elif j > 0 and current_score == score_matrix[i][j-1] + gap_penalty:
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            j -= 1
        # Caso de borda: se o traceback chegar a uma das bordas (i=0 ou j=0)
        else:
            if i > 0: # Se ainda há caracteres em seq2
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = seq2[i-1] + aligned_seq2
                i -= 1
            else: # Se ainda há caracteres em seq1
                aligned_seq1 = seq1[j-1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
                j -= 1
            
    return aligned_seq1, aligned_seq2

if __name__ == "__main__":
    S1 = "GATTACACGA"
    S2 = "GCATGC"

    # O mesmo sistema de pontuação
    MATCH = 1
    MISMATCH = -1
    GAP = -2

    # Gerar a matriz de programação dinâmica usando o código de needleman_wunsch.
    # A função retorna a matriz, as sequências alinhadas e a pontuação.
    matrix, _, _, score = needleman_wunsch(S1, S2, MATCH, MISMATCH, GAP)

    # Obter o alinhamento a partir da matriz gerada.
    aligned_s1_from_matrix, aligned_s2_from_matrix = traceback_from_matrix(S1, S2, matrix, MATCH, MISMATCH, GAP)

    # Exibe os resultados
    print("--- Alinhamento Global a partir da Matriz ---")
    print(f"Sequencia 1 (S1): {S1}")
    print(f"Sequencia 2 (S2): {S2}\n")
    
    print(f"Pontuacao Final: {score}")
    print(f"S1 Alinhada: {aligned_s1_from_matrix}")
    print(f"S2 Alinhada: {aligned_s2_from_matrix}")
