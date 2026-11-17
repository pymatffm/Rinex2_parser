import tempfile
import unittest
from pathlib import Path

from main import (
    open_observation_file,
    read_observation_dataframe,
    summarize_observation_file,
)


def header_line(content, label):
    return "{:<60}{}\n".format(content, label)


def make_rinex_file(contents):
    temporary_file = tempfile.NamedTemporaryFile(
        mode="w", encoding="ascii", suffix=".obs", delete=False
    )
    temporary_file.write(contents)
    temporary_file.close()
    return Path(temporary_file.name)


class RinexObservationParserTests(unittest.TestCase):
    def setUp(self):
        version_line = "{:<60}{}\n".format(
            "     2.11           OBSERVATION DATA    G (GPS)",
            "RINEX VERSION / TYPE",
        )
        observation_types = "{:6d}{}{}\n".format(
            2, "    L1", "    L2"
        )
        self.file_path = make_rinex_file(
            version_line
            + header_line("TEST", "MARKER NAME")
            + header_line(observation_types.rstrip("\n"), "# / TYPES OF OBSERV")
            + header_line("", "END OF HEADER")
            + " 13 11 11  8  0  0.0000000  0  1G01\n"
            + "{:14.3f}{}{}{:14.3f}{}{}\n".format(123.5, 1, 7, 45.25, 0, 6)
        )

    def tearDown(self):
        self.file_path.unlink()

    def test_reads_header_and_observation_fields(self):
        with open_observation_file(self.file_path) as (header, epochs):
            parsed_epochs = list(epochs)

        self.assertEqual(header.version, 2.11)
        self.assertEqual(header.marker_name, "TEST")
        self.assertEqual(header.observation_types, ("L1", "L2"))
        self.assertEqual(len(parsed_epochs), 1)

        epoch = parsed_epochs[0]
        self.assertEqual(epoch.timestamp.isoformat(), "2013-11-11T08:00:00")
        self.assertEqual(epoch.observations["G01"][0].value, 123.5)
        self.assertEqual(epoch.observations["G01"][0].loss_of_lock_indicator, 1)
        self.assertEqual(epoch.observations["G01"][1].signal_strength, 6)

    def test_summary_counts_parsed_values(self):
        summary = summarize_observation_file(self.file_path)

        self.assertIn("Observation epochs: 1", summary)
        self.assertIn("Satellite-epochs: 1", summary)
        self.assertIn("Valid observation values: 2", summary)
        self.assertIn("Satellites: G01", summary)

    def test_dataframe_tracks_each_satellite_at_each_epoch(self):
        header, dataframe = read_observation_dataframe(self.file_path)

        self.assertEqual(len(dataframe), 1)
        self.assertEqual(dataframe.index.names, ["epoch", "satellite"])
        self.assertEqual(dataframe.index[0][1], "G01")
        self.assertEqual(
            dataframe[("L1", "value")].iloc[0], 123.5
        )
        self.assertEqual(dataframe[("L1", "lli")].iloc[0], 1)
        self.assertEqual(dataframe[("L2", "signal_strength")].iloc[0], 6)
        self.assertEqual(
            tuple(dataframe.columns.levels[0]),
            header.observation_types,
        )

    def test_skips_rinex_file_splice_comments(self):
        contents = self.file_path.read_text(encoding="ascii")
        splice = (
            " " * 28 + "4  6\n"
            + header_line("RINEX FILE SPLICE", "COMMENT")
            + header_line("generated comment", "COMMENT")
            + " 13 11 11  8  0  1.0000000  0  1G01\n"
            + "{:14.3f}{}{}{:14.3f}{}{}\n".format(100.0, 0, 5, 50.0, 0, 5)
        )
        spliced_file = make_rinex_file(contents + splice)
        try:
            with open_observation_file(spliced_file) as (header, epochs):
                parsed_epochs = list(epochs)
        finally:
            spliced_file.unlink()

        self.assertEqual(len(parsed_epochs), 2)
        self.assertEqual(parsed_epochs[-1].timestamp.second, 1)


if __name__ == "__main__":
    unittest.main()
