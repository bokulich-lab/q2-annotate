from rachis.plugin.testing import TestPluginBase


class TestUsageExample(TestPluginBase):
    package = "q2_annotate.tests"

    def test_usage_examples(self):
        self.execute_examples()

    def test_tfa_demo_has_scientific_and_partially_resolved_taxonomies(self):
        """Keep valid names and deliberate genus, family, and unclassified cases."""
        import pandas as pd
        from q2_annotate.examples import factory_taxonomy

        species = {
            "T001": "Escherichia coli",
            "T002": "Klebsiella pneumoniae",
            "T003": "Enterococcus faecium",
            "T004": "Bacteroides fragilis",
            "T005": "Phocaeicola vulgatus",
            "T006": "Faecalibacterium prausnitzii",
            "T007": "Roseburia intestinalis",
            "T008": "Akkermansia muciniphila",
            "T009": "Bifidobacterium longum",
            "T011": "Clostridioides difficile",
        }
        expected = {
            identifier: {
                "unclassified": False,
                "ranks": list("dpcofgs"),
                "all_named": True,
                "genus": name.split()[0],
                "species": name,
            }
            for identifier, name in species.items()
        }
        expected.update(
            {
                "T010": {
                    "unclassified": False,
                    "ranks": list("dpcofg"),
                    "all_named": True,
                    "genus": "Segatella",
                    "species": None,
                },
                "T012": {
                    "unclassified": False,
                    "ranks": list("dpcof"),
                    "all_named": True,
                    "genus": None,
                    "species": None,
                },
                "T013": {
                    "unclassified": True,
                    "ranks": [],
                    "all_named": False,
                    "genus": None,
                    "species": None,
                },
            }
        )
        taxonomy = factory_taxonomy().view(pd.DataFrame)["Taxon"]
        observed = {}
        for identifier, label in taxonomy.items():
            unclassified = label == "d__Unclassified"
            ranks = (
                []
                if unclassified
                else [rank.strip().split("__", 1) for rank in label.split(";")]
            )
            names = dict(ranks)
            observed[identifier] = {
                "unclassified": unclassified,
                "ranks": [rank for rank, _ in ranks],
                "all_named": bool(ranks) and all(name for _, name in ranks),
                "genus": names.get("g"),
                "species": names.get("s"),
            }
        self.assertEqual(observed, expected)
