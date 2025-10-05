import numpy as np

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    """
    Implementa o algoritmo Needleman-Wunsch para alinhamento global de sequências.

    Args:
        seq1 (str): A primeira sequência.
        seq2 (str): A segunda sequência.
        match_score (int): Pontuação para um match.
        mismatch_score (int): Pontuação para um mismatch.
        gap_penalty (int): Pontuação para um gap (lacuna).

    Returns:
        tuple: Uma tupla contendo a matriz de pontuação (numpy.ndarray),
               a primeira sequência alinhada (str),
               a segunda sequência alinhada (str),
               e a pontuação final do alinhamento (int).
    """
    # Obtém os tamanhos das sequências
    n = len(seq1)
    m = len(seq2)

    # Inicializa a matriz de pontuação com zeros.
    # As dimensões são (m+1) x (n+1) para acomodar os gaps iniciais.
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)

    # --- 1. Inicialização da Matriz ---
    # A primeira linha e a primeira coluna são preenchidas com as penalidades de gap.
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i-1][0] + gap_penalty
    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j-1] + gap_penalty

    # --- 2. Preenchimento da Matriz ---
    # Itera sobre cada célula da matriz para calcular a pontuação ótima.
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calcula a pontuação para um match ou mismatch (vindo da diagonal)
            char1 = seq1[j-1]
            char2 = seq2[i-1]
            match_value = match_score if char1 == char2 else mismatch_score
            diagonal_score = score_matrix[i-1][j-1] + match_value

            # Calcula a pontuação para um gap na seq1 (vindo de cima)
            up_score = score_matrix[i-1][j] + gap_penalty

            # Calcula a pontuação para um gap na seq2 (vindo da esquerda)
            left_score = score_matrix[i][j-1] + gap_penalty

            # A pontuação da célula atual é o máximo das três possibilidades
            score_matrix[i][j] = max(diagonal_score, up_score, left_score)

    # A pontuação final do alinhamento está na célula inferior direita
    alignment_score = score_matrix[m][n]

    # --- 3. Traceback (Rastreamento do Caminho) ---
    # Começa do canto inferior direito da matriz para construir as sequências alinhadas.
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = m, n

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        # Verifica se o caminho veio da diagonal
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
        # O caminho veio da esquerda (gap em seq2)
        else:
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            j -= 1
            
    return score_matrix, aligned_seq1, aligned_seq2, alignment_score

# --- Execução Principal ---
if __name__ == "__main__":
    # Sequências de DNA de exemplo
    S1 = "GATTACACGA"
    S2 = "GCATGC"

    # Sistema de pontuação
    MATCH = 1
    MISMATCH = -1
    GAP = -2

    # Executa o algoritmo
    matrix, aligned_s1, aligned_s2, score = needleman_wunsch(S1, S2, MATCH, MISMATCH, GAP)

    # Exibe os resultados
    print("--- Algoritmo Needleman-Wunsch ---")
    print(f"Sequencia 1 (S1): {S1}")
    print(f"Sequencia 2 (S2): {S2}\n")

    print("Sistema de Pontuacao:")
    print(f"  - Match: {MATCH}")
    print(f"  - Mismatch: {MISMATCH}")
    print(f"  - Gap: {GAP}\n")

    print("--- Matriz de Programacao Dinamica Gerada ---")
    # Adiciona cabeçalhos para melhor visualização
    header_s1 = "      " + "   ".join(S1)
    print(header_s1)
    for i in range(matrix.shape[0]):
        row_header = " " if i == 0 else S2[i-1]
        print(f"{row_header} {matrix[i,:]}")
    
    print("\n--- Resultado do Alinhamento Global ---")
    print(f"Pontuacao Final: {score}")
    print(f"S1 Alinhada: {aligned_s1}")
    print(f"S2 Alinhada: {aligned_s2}")
