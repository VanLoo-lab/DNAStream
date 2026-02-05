from pathlib import Path
from dnastream._builtin_schemas import REGISTRY_SCHEMAS, PROVENANCE_SCHEMAS

out = []
out.append("# Built-in schemas\n")
out.append("These tables are generated from the current `dnastream` code.\n")

out.append("## Registries\n")
for name, schema in REGISTRY_SCHEMAS.items():
    out.append(f"### `{name}`\n")
    out.append(schema.to_markdown_table())
    out.append("")

out.append("## Provenance\n")
for name, schema in PROVENANCE_SCHEMAS.items():
    out.append(f"### `{name}`\n")
    out.append(schema.to_markdown_table())
    out.append("")

Path("docs/reference/schemas.md").write_text("\n".join(out), encoding="utf-8")
print("Wrote docs/reference/schemas.md")
