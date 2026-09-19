"""Tournaments: a simple directed graph where every pair of distinct
vertices has exactly one directed edge between them (the result of a
round-robin contest). Demonstrates: transitivity testing, strongly-
connected-component decomposition purely from in-degrees, a Hamiltonian
path (Redei's theorem), a Hamiltonian cycle for strongly connected
tournaments (Camion's theorem), and a 2-king (Landau's theorem).
"""
from __future__ import annotations

import random


def random_tournament(n: int, rng: random.Random) -> list[list[bool]]:
    """beats[i][j] = True iff i -> j (i beat j)."""
    beats = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                beats[i][j] = True
            else:
                beats[j][i] = True
    return beats


def is_transitive(beats: list[list[bool]], n: int) -> bool:
    """Equivalent to "contains no 3-cycle"."""
    for a in range(n):
        for b in range(n):
            if not beats[a][b]:
                continue
            for c in range(n):
                if beats[b][c] and beats[c][a]:
                    return False
    return True


def in_degrees(beats: list[list[bool]], n: int) -> list[int]:
    return [sum(1 for u in range(n) if beats[u][v]) for v in range(n)]


def scc_by_reachability(beats: list[list[bool]], n: int) -> list[list[int]]:
    """Reference SCC decomposition via mutual reachability (O(n^3)
    transitive closure), used only to cross-check the in-degree method."""
    reach = [[False] * n for _ in range(n)]
    for i in range(n):
        reach[i][i] = True
        for j in range(n):
            if beats[i][j]:
                reach[i][j] = True
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]:
                        reach[i][j] = True

    assigned = [False] * n
    comps = []
    for v in range(n):
        if assigned[v]:
            continue
        comp = [u for u in range(n) if reach[v][u] and reach[u][v]]
        for u in comp:
            assigned[u] = True
        comps.append(comp)
    return comps


def scc_by_indegree(beats: list[list[bool]], n: int) -> list[list[int]]:
    """SCC decomposition using ONLY in-degrees: sort vertices ascending by
    in-degree, then cut wherever the running sum of in-degrees among the
    first k equals C(k,2) -- the signature of "no edges enter from
    outside".
    """
    indeg = in_degrees(beats, n)
    order = sorted(range(n), key=lambda v: indeg[v])

    comps = []
    start = 0
    running = 0
    for k in range(1, n + 1):
        running += indeg[order[k - 1]]
        target = k * (k - 1) // 2
        if running == target:
            comps.append(order[start:k])
            start = k
    return comps


def hamiltonian_path(beats: list[list[bool]], n: int) -> list[int]:
    """Redei's theorem: build a Hamiltonian path by inserting vertices one
    at a time into the correct position of the path built so far.
    """
    path: list[int] = []
    for v in range(n):
        if not path:
            path.append(v)
            continue
        if beats[v][path[0]]:
            path.insert(0, v)
        elif beats[path[-1]][v]:
            path.append(v)
        else:
            idx = 0
            for i in range(len(path) - 1):
                if beats[path[i]][v] and beats[v][path[i + 1]]:
                    idx = i + 1
                    break
            path.insert(idx, v)
    return path


def is_valid_hamiltonian_path(beats: list[list[bool]], path: list[int], n: int) -> bool:
    if len(path) != n or len(set(path)) != n:
        return False
    return all(beats[path[i]][path[i + 1]] for i in range(n - 1))


def hamiltonian_cycle(beats: list[list[bool]], n: int, path: list[int]) -> list[int]:
    """Camion's theorem: for a strongly connected tournament with n >= 3,
    extend a Hamiltonian path's growing prefix cycle one step at a time
    until it covers everything.
    """
    order = path[:]
    k = 0
    for j in range(n - 1, 0, -1):
        if beats[order[j]][order[0]]:
            k = j + 1
            break

    while k < n:
        next_v = order[k]  # this is v_{k+1}
        if beats[next_v][order[0]]:
            k += 1  # Case A: extends directly.
            continue

        found_p = None
        for p0 in range(k):
            if beats[next_v][order[p0]]:
                found_p = p0
                break

        if found_p is not None:
            p0 = found_p  # p0 >= 1 always here
            order = order[p0:k] + order[0:p0] + [next_v] + order[k + 1 : n]
            k += 1
        else:
            # Case C: every current cycle vertex beats next_v. Use strong
            # connectivity to find the first later vertex with an edge
            # back into the cycle.
            found = None
            for q0 in range(k + 1, n):
                for p0 in range(k):
                    if beats[order[q0]][order[p0]]:
                        found = (q0, p0)
                        break
                if found:
                    break
            assert found is not None, "strongly connected tournament guarantees this exists"
            q0, p0 = found
            if p0 == 0:
                k = q0 + 1  # order[q0] -> order[0] directly: no reordering needed
            else:
                order = order[p0:k] + order[0:p0] + order[k : q0 + 1] + order[q0 + 1 : n]
                k = q0 + 1
    return order


def is_valid_hamiltonian_cycle(beats: list[list[bool]], cycle: list[int], n: int) -> bool:
    if not is_valid_hamiltonian_path(beats, cycle, n):
        return False
    return beats[cycle[n - 1]][cycle[0]]


def find_2_king(beats: list[list[bool]], n: int) -> int:
    """Landau's theorem: the maximum out-degree vertex reaches every other
    vertex within distance 2."""
    return max(range(n), key=lambda v: sum(1 for u in range(n) if beats[v][u]))


def is_2_king(beats: list[list[bool]], v: int, n: int) -> bool:
    return all(
        x == v or beats[v][x] or any(beats[v][y] and beats[y][x] for y in range(n))
        for x in range(n)
    )


def main() -> None:
    rng = random.Random(12345)

    # --- Transitivity, SCC decomposition, in-degree cross-check --------
    scc_ok = True
    for trial in range(200):
        n = 2 + (trial % 6)
        beats = random_tournament(n, rng)
        ref_set = sorted(sorted(c) for c in scc_by_reachability(beats, n))
        deg_set = sorted(sorted(c) for c in scc_by_indegree(beats, n))
        if ref_set != deg_set:
            scc_ok = False
            print(f"SCC MISMATCH n={n}: {ref_set} vs {deg_set}")
    print(f"in-degree-only SCC decomposition matches reachability-based SCC: {'PASS' if scc_ok else 'FAIL'}")

    # --- Hamiltonian path, on many random tournaments -------------------
    path_ok = True
    for trial in range(200):
        n = 1 + (trial % 8)
        beats = random_tournament(n, rng)
        path = hamiltonian_path(beats, n)
        if not is_valid_hamiltonian_path(beats, path, n):
            path_ok = False
            print(f"HAM PATH FAIL n={n} beats={beats} path={path}")
    print(f"Redei Hamiltonian path valid on every trial: {'PASS' if path_ok else 'FAIL'}")

    # --- Hamiltonian cycle, on strongly connected tournaments ------------
    cycle_ok = True
    cycle_trials = 0
    trial = 0
    while cycle_trials < 200:
        trial += 1
        n = 3 + (trial % 6)
        beats = random_tournament(n, rng)
        if len(scc_by_reachability(beats, n)) != 1:
            continue
        cycle_trials += 1
        path = hamiltonian_path(beats, n)
        cycle = hamiltonian_cycle(beats, n, path)
        if not is_valid_hamiltonian_cycle(beats, cycle, n):
            cycle_ok = False
            print(f"HAM CYCLE FAIL n={n} beats={beats} cycle={cycle}")
    print(f"Camion Hamiltonian cycle valid on every strongly-connected trial: {'PASS' if cycle_ok else 'FAIL'}")

    # --- 2-king ------------------------------------------------------------
    king_ok = True
    for trial in range(200):
        n = 1 + (trial % 8)
        beats = random_tournament(n, rng)
        king = find_2_king(beats, n)
        if not is_2_king(beats, king, n):
            king_ok = False
            print(f"2-KING FAIL n={n} beats={beats} king={king}")
    print(f"max-out-degree vertex is always a 2-king: {'PASS' if king_ok else 'FAIL'}")

    # --- A concrete worked example: rock-paper-scissors -------------------
    # 0=Rock 1=Scissors 2=Paper (Rock beats Scissors, Scissors beats Paper,
    # Paper beats Rock): a 3-cycle, so NOT transitive.
    n = 3
    beats = [[False] * n for _ in range(n)]
    beats[0][1] = True  # Rock beats Scissors
    beats[1][2] = True  # Scissors beats Paper
    beats[2][0] = True  # Paper beats Rock
    print(f"\nrock-paper-scissors tournament transitive? {is_transitive(beats, n)}")
    print(f"its SCCs (in-degree method): {scc_by_indegree(beats, n)}")


if __name__ == "__main__":
    main()
