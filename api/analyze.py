
from http.server import BaseHTTPRequestHandler
import json

ALLOWED = set("ATGC")


def validate(seq, name="sequence"):
    if not isinstance(seq, str):
        raise ValueError(f"{name} must be text.")

    seq = "".join(seq.upper().split())

    if not seq:
        raise ValueError(f"{name} cannot be empty.")

    if set(seq) - ALLOWED:
        raise ValueError("Only A, T, G and C are allowed.")

    if len(seq) > 10000:
        raise ValueError("Maximum sequence length is 10000.")

    return seq


def align_sequences(a, b):
    match, mismatch, gap = 2, -1, -2
    n, m = len(a), len(b)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * gap
        trace[i][0] = "U"

    for j in range(1, m + 1):
        dp[0][j] = j * gap
        trace[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = dp[i-1][j-1] + (
                match if a[i-1] == b[j-1] else mismatch
            )
            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap

            best = max(diagonal, up, left)
            dp[i][j] = best

            if best == diagonal:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    x, y = [], []
    i, j = n, m

    while i > 0 or j > 0:
        direction = trace[i][j]

        if direction == "D":
            x.append(a[i-1])
            y.append(b[j-1])
            i -= 1
            j -= 1

        elif direction == "U":
            x.append(a[i-1])
            y.append("-")
            i -= 1

        else:
            x.append("-")
            y.append(b[j-1])
            j -= 1

    return "".join(reversed(x)), "".join(reversed(y))


def compare(reference, sample):
    a, b = align_sequences(reference, sample)

    mutations = []
    ref_pos = 0
    sample_pos = 0

    matches = 0
    substitutions = 0
    insertions = 0
    deletions = 0

    for x, y in zip(a, b):

        if x != "-":
            ref_pos += 1

        if y != "-":
            sample_pos += 1

        if x == y:
            matches += 1

        elif x != "-" and y != "-":
            substitutions += 1
            mutations.append({
                "type": "Substitution",
                "reference": x,
                "sample": y,
                "position": ref_pos
            })

        elif x == "-":
            insertions += 1
            mutations.append({
                "type": "Insertion",
                "reference": "",
                "sample": y,
                "position": sample_pos
            })

        else:
            deletions += 1
            mutations.append({
                "type": "Deletion",
                "reference": x,
                "sample": "",
                "position": ref_pos
            })

    gc = (
        (sample.count("G") + sample.count("C"))
        / len(sample) * 100
    )

    return {
        "summary": {
            "reference_length": len(reference),
            "sample_length": len(sample),
            "matches": matches,
            "substitutions": substitutions,
            "insertions": insertions,
            "deletions": deletions,
            "gc_content_percent": round(gc, 2)
        },
        "mutations": mutations,
        "alignment": {
            "reference": a,
            "sample": b
        }
    }


def find_pattern(sequence, pattern):
    positions = []
    start = 0

    while True:
        index = sequence.find(pattern, start)

        if index == -1:
            break

        positions.append(index + 1)
        start = index + 1

    return {
        "count": len(positions),
        "positions": positions,
        "overlapping": True
    }


def edit_sequence(sequence, operation, position, base):
    if not isinstance(position, int) or position < 1:
        raise ValueError("Position must be a positive integer.")

    if operation not in ("replace", "insert", "delete"):
        raise ValueError("Invalid operation.")

    if operation == "insert":
        if position > len(sequence) + 1:
            raise ValueError("Position is outside the sequence.")

    elif position > len(sequence):
        raise ValueError("Position is outside the sequence.")

    if operation in ("replace", "insert"):
        if base not in ALLOWED:
            raise ValueError("Invalid nucleotide.")

    i = position - 1

    if operation == "replace":
        modified = sequence[:i] + base + sequence[i+1:]
        message = f"Position {position} replaced with {base}."

    elif operation == "insert":
        modified = sequence[:i] + base + sequence[i:]
        message = f"{base} inserted before position {position}."

    else:
        modified = sequence[:i] + sequence[i+1:]
        message = f"Position {position} deleted."

    return {
        "original": sequence,
        "modified": modified,
        "message": message
    }


class handler(BaseHTTPRequestHandler):

    def send_json(self, status, data):
        body = json.dumps(data).encode()

        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_json(200, {"success": True})

    def do_POST(self):
        try:
            length = int(self.headers.get("Content-Length", 0))

            if length > 100000:
                return self.send_json(413, {
                    "error": "Request too large."
                })

            data = json.loads(self.rfile.read(length))
            action = data.get("action")

            if action == "compare":
                reference = validate(
                    data.get("reference"), "Reference sequence"
                )
                sample = validate(
                    data.get("sample"), "Sample sequence"
                )

                self.send_json(
                    200, compare(reference, sample)
                )

            elif action == "find":
                sequence = validate(data.get("sequence"))
                pattern = validate(data.get("pattern"), "Pattern")

                self.send_json(
                    200, find_pattern(sequence, pattern)
                )

            elif action == "edit":
                sequence = validate(data.get("sequence"))

                result = edit_sequence(
                    sequence,
                    data.get("type"),
                    data.get("position"),
                    data.get("base", "").upper()
                )

                self.send_json(200, result)

            else:
                self.send_json(400, {
                    "error": "Unknown operation."
                })

        except (ValueError, json.JSONDecodeError) as e:
            self.send_json(400, {"error": str(e)})

        except Exception:
            self.send_json(500, {
                "error": "Internal server error."
            })
