from pathlib import Path

datasets = {
    file.name.replace('.relative_abundance.tsv', ''): file
    for file in Path(__file__).parent.rglob('*.relative_abundance.tsv')
}
