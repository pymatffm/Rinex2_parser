import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rinex_2_tool import (
    calculate_tec_dataframe,
    main,
    open_observation_file,
    read_observation_dataframe,
    write_azimuth_elevation_files,
)
from rinex_visualization import (
    plot_pierce_point_map,
    plot_vtec_by_satellite,
)
from tec_calculation import (
    calculate_tec_dataframe as calculate_normalized_tec_dataframe,
    normalize_long_observations,
    observation_time_to_utc,
    stec_from_carrier_phase,
    stec_from_pseudoranges,
    thin_shell_mapping_factor,
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

    def test_code_tec_and_mapping_factor(self):
        frequency_1 = 1575.42e6
        frequency_2 = 1227.60e6
        self.assertEqual(
            stec_from_pseudoranges(
                20_000_000.0, 20_000_000.0, frequency_1, frequency_2
            ),
            0.0,
        )
        self.assertGreater(
            stec_from_pseudoranges(
                20_000_000.0, 20_000_010.0, frequency_1, frequency_2
            ),
            0.0,
        )
        self.assertGreater(
            stec_from_carrier_phase(100.0, 70.0, frequency_1, frequency_2), 0.0
        )
        self.assertGreater(thin_shell_mapping_factor(90.0, 350.0), 0.99)
        self.assertLess(thin_shell_mapping_factor(5.0, 350.0), 0.4)

    def test_gps_observation_time_is_converted_to_utc(self):
        gps_time = datetime(2013, 11, 11, 8)
        self.assertEqual(
            observation_time_to_utc(gps_time, "GPS"),
            datetime(2013, 11, 11, 7, 59, 44),
        )
        self.assertEqual(
            observation_time_to_utc(gps_time, "BDT"),
            datetime(2013, 11, 11, 7, 59, 58),
        )
        self.assertEqual(
            observation_time_to_utc(gps_time, "TAI"),
            datetime(2013, 11, 11, 7, 59, 25),
        )
        self.assertEqual(observation_time_to_utc(gps_time, "UTC"), gps_time)

    def test_normalized_tec_api_accepts_other_signal_frequencies(self):
        first_epoch = datetime(2013, 11, 11, 8)
        observations = pd.DataFrame([
            {
                "epoch": first_epoch + pd.Timedelta(seconds=second),
                "satellite": "E01",
                "code_1_m": 22_000_000.0,
                "code_2_m": 22_000_005.0 + second * 0.1,
                "phase_1_cycles": 100_000.0 + second,
                "phase_2_cycles": 80_000.0 + second,
                "frequency_1_hz": 1575.42e6,
                "frequency_2_hz": 1176.45e6,
                "lli_1": 0,
                "lli_2": 0,
                "code_pair": "C1C/C5Q",
            }
            for second in (0, 1, 2)
        ])
        elevations = pd.DataFrame([
            {
                "epoch": first_epoch + pd.Timedelta(seconds=second),
                "satellite": "E01",
                "elevation_deg": 45.0,
            }
            for second in (0, 1, 2)
        ])

        result = calculate_normalized_tec_dataframe(
            observations, elevations, time_system="GST"
        )

        self.assertEqual(len(result), 3)
        self.assertEqual(set(result["satellite"]), {"E01"})
        self.assertEqual(result["code_pair"].iloc[0], "C1C/C5Q")
        self.assertEqual(result["epoch"].iloc[0], datetime(2013, 11, 11, 7, 59, 44))

    def test_rinex34_long_rows_can_be_normalized_for_tec(self):
        epoch = datetime(2024, 3, 14)
        long_observations = pd.DataFrame([
            {
                "time": epoch,
                "sat": "E11",
                "code": code,
                "value": value,
                "lli": 0,
            }
            for code, value in (
                ("C1C", 22_000_000.0),
                ("C5Q", 22_000_004.0),
                ("L1C", 100_000.0),
                ("L5Q", 80_000.0),
            )
        ])

        normalized = normalize_long_observations(
            long_observations,
            {
                "E": (
                    "C1C",
                    "C5Q",
                    "L1C",
                    "L5Q",
                    1575.42e6,
                    1176.45e6,
                )
            },
        )
        elevations = pd.DataFrame(
            [{"epoch": epoch, "satellite": "E11", "elevation_deg": 45.0}]
        )
        result = calculate_normalized_tec_dataframe(
            normalized, elevations, time_system="GPS"
        )

        self.assertEqual(result["satellite"].tolist(), ["E11"])
        self.assertEqual(result["code_pair"].tolist(), ["C1C/C5Q"])

    def test_pierce_point_map_reads_compact_angle_files(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            azimuth_path = directory / "sample.azi"
            elevation_path = directory / "sample.ele"
            header = (
                "COMPACT2\n"
                "T_SAMP    15.0000\n"
                "START_TIME_MJD  56609.000000000\n"
                "1 G01\n"
                "  90.000\n"
                "-1\n"
                "  91.000\n"
            )
            azimuth_path.write_text(header, encoding="ascii")
            elevation_path.write_text(
                header.replace("  90.000", "  45.000").replace(
                    "  91.000", "  46.000"
                ),
                encoding="ascii",
            )

            output_path = plot_pierce_point_map(
                azimuth_path,
                elevation_path,
                receiver_latitude=-60.5,
                receiver_longitude=-50.0,
            )

            self.assertTrue(output_path.is_file())
            self.assertGreater(output_path.stat().st_size, 0)

    def test_tec_dataframe_uses_phase_smoothing_and_creates_plot(self):
        repository_root = Path(__file__).resolve().parents[1]
        observation_path = repository_root / "Data_2" / "kepa3170.13o"
        elevation_path = repository_root / "Data_2" / "kepa3170.ele"
        dataframe = calculate_tec_dataframe(observation_path, elevation_path)

        self.assertGreater(len(dataframe), 0)
        self.assertEqual(
            dataframe["epoch"].iloc[0], datetime(2013, 11, 13, 7, 59, 44)
        )
        self.assertIn("C1/C2", set(dataframe["code_pair"]))
        self.assertTrue(
            dataframe["phase_smoothed_vtec_tecu_uncalibrated"].notna().all()
        )
        self.assertFalse(
            dataframe["phase_smoothed_vtec_tecu_uncalibrated"].equals(
                dataframe["vtec_tecu_uncalibrated"]
            )
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            plots = plot_vtec_by_satellite(dataframe, Path(temporary_directory))
            self.assertEqual(len(plots), dataframe["satellite"].nunique())
            self.assertTrue(all(path.is_file() and path.stat().st_size > 0 for path in plots))

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

    def test_skips_embedded_header_metadata_between_epochs(self):
        contents = self.file_path.read_text(encoding="ascii")
        spliced_file = make_rinex_file(
            contents
            + header_line("updated receiver position", "COMMENT")
            + header_line("1.0 2.0 3.0", "APPROX POSITION XYZ")
            + " 13 11 11  8  0  1.0000000  0  1G01\n"
            + "{:14.3f}{}{}{:14.3f}{}{}\n".format(100.0, 0, 5, 50.0, 0, 5)
        )
        try:
            with open_observation_file(spliced_file) as (header, epochs):
                parsed_epochs = list(epochs)
        finally:
            spliced_file.unlink()

        self.assertEqual(len(parsed_epochs), 2)
        self.assertEqual(parsed_epochs[-1].timestamp.second, 1)

    def test_cli_writes_csv_into_output_directory_by_default(self):
        output_dir = Path(tempfile.mkdtemp())
        try:
            with patch("rinex_2_tool.DATA_2_OUTPUT_DIR", output_dir):
                result = main([str(self.file_path)])

            output_file = output_dir / (self.file_path.name + ".csv")
            self.assertEqual(result, 0)
            self.assertTrue(output_file.is_file())
            self.assertIn("G01", output_file.read_text(encoding="utf-8"))
        finally:
            for output_file in output_dir.iterdir():
                output_file.unlink()
            output_dir.rmdir()

    def test_navigation_cli_writes_csv_into_output_directory(self):
        fields = [1.0] * 29
        first_line = "{:2d}{:3d}{:3d}{:3d}{:3d}{:3d}{:5.1f}".format(
            1, 13, 11, 11, 8, 0, 0.0
        )
        first_line += "".join(
            "{:19.12E}".format(value).replace("E", "D")
            for value in fields[:3]
        )
        continuation_lines = []
        for offset in range(3, len(fields), 4):
            line_fields = fields[offset:offset + 4]
            continuation_lines.append(
                "   "
                + "".join(
                    "{:19.12E}".format(value).replace("E", "D")
                    for value in line_fields
                )
            )
        contents = (
            header_line(
                "     2.11           NAVIGATION DATA    G (GPS)",
                "RINEX VERSION / TYPE",
            )
            + header_line("", "END OF HEADER")
            + first_line
            + "\n"
            + "\n".join(continuation_lines)
            + "\n"
        )
        navigation_file = make_rinex_file(contents)
        output_dir = Path(tempfile.mkdtemp())
        try:
            with patch("rinex_2_tool.DATA_2_OUTPUT_DIR", output_dir):
                result = main(["nav", str(navigation_file)])

            output_file = output_dir / (navigation_file.name + ".csv")
            self.assertEqual(result, 0)
            self.assertTrue(output_file.is_file())
            self.assertIn("G01", output_file.read_text(encoding="utf-8"))
        finally:
            navigation_file.unlink()
            for output_file in output_dir.iterdir():
                output_file.unlink()
            output_dir.rmdir()

    def test_angles_match_teqc_compact_reference(self):
        repository_root = Path(__file__).resolve().parents[1]
        data_dir = repository_root / "Data_2"
        observation_path = data_dir / "kepa3170.13o"
        navigation_path = data_dir / "brdc3170.13n"
        reference_paths = {
            "azi": data_dir / "kepa3170.azi",
            "ele": data_dir / "kepa3170.ele",
        }

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_paths = write_azimuth_elevation_files(
                observation_path,
                navigation_path,
                Path(temporary_directory),
            )

            def read_compact_rows(path):
                lines = path.read_text(encoding="ascii").splitlines()
                satellite_ids = []
                rows = []
                for line_number in range(3, len(lines), 2):
                    marker = lines[line_number].split()
                    if marker[0] != "-1":
                        satellite_ids = marker[1:]
                    rows.append(
                        dict(zip(satellite_ids, map(float, lines[line_number + 1].split())))
                    )
                return rows

            generated_azi = read_compact_rows(output_paths[0])
            generated_ele = read_compact_rows(output_paths[1])
            reference_azi = read_compact_rows(reference_paths["azi"])
            reference_ele = read_compact_rows(reference_paths["ele"])

        self.assertEqual(len(generated_azi), 480)
        self.assertEqual(len(generated_ele), 480)
        azimuth_errors = []
        elevation_errors = []
        for offset, (azi_row, ele_row) in enumerate(
            zip(generated_azi, generated_ele)
        ):
            reference_index = 1920 + offset
            expected_azi = reference_azi[reference_index]
            expected_ele = reference_ele[reference_index]
            for satellite in set(azi_row) ^ set(expected_azi):
                elevation = (
                    ele_row[satellite]
                    if satellite in ele_row
                    else expected_ele[satellite]
                )
                self.assertLessEqual(abs(elevation), 0.02)
            self.assertEqual(set(azi_row), set(ele_row))
            common_satellites = set(azi_row) & set(expected_azi)
            azimuth_errors.extend(
                abs((azi_row[satellite] - expected_azi[satellite] + 180) % 360 - 180)
                for satellite in common_satellites
            )
            elevation_errors.extend(
                abs(ele_row[satellite] - expected_ele[satellite])
                for satellite in common_satellites
            )

        self.assertLessEqual(max(azimuth_errors), 0.07)
        self.assertLessEqual(max(elevation_errors), 0.015)


if __name__ == "__main__":
    unittest.main()
