from ..dnastream import DNAStream


def add_parser(subparsers):
    p = subparsers.add_parser("create", help="Create a DNAStream file")
    p.add_argument(
        "-f", "--file", required=True, help="Path to the DNAStream file to be created."
    )
    p.add_argument(
        "-p",
        "--patient-id",
        required=False,
        help="Patient identifier to store in file (optional).",
    )
    p.add_argument("-v", "--verbose", action="store_true")
    p.set_defaults(func=run)
    return p


def run(args) -> int:
    ds = DNAStream.create(args.file, patient_id=args.patient_id, verbose=args.verbose)
    try:
        return 0
    finally:
        ds.close()
