# Publication bot

Keeps the [Publications](../docs/src/man/publications.md) page of the documentation up to date with
papers that use GeophysicalModelGenerator.jl. Adapted from the publication bot of
[GPU4GEO](https://github.com/PTsolvers/GPU4GEO) (MIT license, see [LICENSE](LICENSE)). Two GitHub
workflows use it:

- **Scan** ([PubScan.yml](../.github/workflows/PubScan.yml)): on the 1st and 15th of each month,
  finds recent papers that use the package and opens a pull request that adds them. A scan is
  skipped while the previous pull request is open.
- **Add a publication** issue form ([PubAdd.yml](../.github/workflows/PubAdd.yml)): repository
  collaborators paste DOIs, and a pull request with the formatted entries is opened.

## How papers are selected

1. [OpenAlex](https://openalex.org) works published within the last `lookback_months` that name the
   package in their full text (OpenAlex's index), cite a package paper (`dois`), or are written
   by the `developers` ([config.toml](config.toml)).
2. Repository copies, peer-review comments and duplicates are dropped; a preprint and its published
   version are merged. DOIs already on the page or listed in `seen_dois.txt` are skipped.
3. Each paper is classified from its title, abstract, open-access PDF and the OpenAlex full-text index:
   - ✅ **added**: the title, abstract or PDF names GeophysicalModelGenerator, whoever the authors are;
   - 🟡 **needs review**: only the OpenAlex full-text index finds the name (it matches word
     stems), or it only cites a package paper (it may be cited as related work); listed in the
     pull request with a ready-to-paste entry;
   - ❌ **no signal**: listed collapsed in the pull request.
4. Preprints on the page whose published version appeared are replaced.

Merging a scan pull request records every evaluated DOI in `seen_dois.txt`, so nothing is proposed
twice. To drop an added entry, delete its line in the pull request before merging; to add a paper
from the review list, paste its entry under the right year. To re-evaluate a paper, delete its DOI
from `seen_dois.txt`.

## Running locally

```sh
julia --project=_pubbot -e 'using Pkg; Pkg.instantiate()'
julia --project=_pubbot _pubbot/weekly.jl --dry-run            # print the report only
julia --project=_pubbot _pubbot/add_doi.jl 10.1038/s41467-024-55443-z
julia --project=_pubbot -e 'using Pkg; Pkg.test()'
```

PDF text is cached in `_pubbot/.cache` (git-ignored).

## OpenAlex API key

OpenAlex answers queries with more than one page of results (50 works) only with an API key, and
without one the runner shares a small anonymous daily budget with everyone on its IP address. Store a
free key from <https://openalex.org/settings/api> as the repository secret `OPENALEX_API_KEY`;
locally, export the same variable.
