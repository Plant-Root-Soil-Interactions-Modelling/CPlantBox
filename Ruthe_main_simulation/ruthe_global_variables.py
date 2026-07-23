"""Configuration for Ruthe growth periods using dataclasses.

This module provides a dataclass-based configuration (no validation) and
multiple instances for different growth periods.

Usage examples:
        from ruthe_global_variables import RUTHE_1995_96

        # Backward-compatible module constants still exist:
        # from ruthe_global_variables import NUMBER_OF_ROWS
"""

from dataclasses import dataclass
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_RUTHE_RAW = BASE_DIR / "data" / "data_raw"
DATA_RUTHE_PROCESSED = BASE_DIR / "data" / "data_processed"
DEFAULT_OUTPUT_DIR = BASE_DIR / "results"


@dataclass
class RutheConfig:
    """Dataclass holding configuration for all Ruthe growth periods."""

    # Parameters of experimental field
    NUMBER_OF_ROWS_SOIL_CORE: int = (
        6  # number of plant rows in the field considered in soil core simulation
    )
    NUMBER_OF_PLANTS_PER_ROW_SOIL_CORE: int = (
        4  # number of plants per row in the field considered in soil core simulation
    )
    NUMBER_OF_ROWS_MINIRHIZOTRON: int = (
        3  # number of plant rows in the field considered in minirhizotron simulation
    )
    NUMBER_OF_PLANTS_PER_ROW_MINIRHIZOTRON: int = 60  # number of plants per row in the field considered in minirhizotron simulation
    INTER_PLANT_DISTANCE: float = 2.27  # cm, estimated distance between root systems
    INTER_ROW_DISTANCE: float = 12.5  # cm, distance between plant rows in the field
    SOIL_SPACE_PARAMS: tuple[int, int, int, bool] = (500, 500, 500, True)
    SOWING_DEPTH: float = -3.0

    # Soil temperature range and optimum for scaling the root elongation
    TEMP_MIN: float = 2.0  # minimum temperature required for growth
    TEMP_MAX: float = 25.0  # maximum temperature for growth
    TEMP_OPT: float = 16.3  # optimal temperature for growth

    # Parameters for soil cores taken
    SOIL_CORE_RADIUS: float = 4.0  # cm
    SOIL_CORE_DEPTH: float = 120.0  # cm
    NUMBER_OF_LAYERS: int = 8  # number of soil layers in the soil core; 15cm thick each
    NUMBER_OF_SOIL_CORE_OBSERVATIONS: int = (
        6  # number of soil core samples taken per depth
    )
    MEASUREMENT_DAY_SOIL_CORE: str = (None)

    # Parameters needed for minirhizotron simulation
    TUBE_EXTERNAL_DIAMETER: float = (
        5  # cm  --> viewing depth = (external_diameter - internal_diameter) / 2
    )
    TUBE_INTERNAL_DIAMETER: float = 4.6  # cm
    TUBE_OPENING_DIAMETER: float = 3.0  # cm
    TUBE_OPENING_HEIGHT: float = 5.0  # cm
    INTER_IMAGE_DISTANCE: float = -5.77  # cm (negative value indicates depth direction)
    TUBE_ANGLE: float = 30.0  # degrees from vertical
    TUBE_LENGTH: float = 180.0  # cm
    NUMBER_OF_LAYERS_MINIRHIZOTRON: int = (
        18  # number of soil layers in the minirhizotron; 15cm thick each
    )
    NUMBER_OF_MINIRHIZO_OBSERVATIONS: int = (
        24  # number of minirhizotron samples taken per depth
    )
    BONITUR_SCORING_BASE: float = (
        0.5  # base for scoring minirhizotron root observations
    )
    DAYS_TO_ANALYZE_MINIRHIZOTRON: list[int] = (
        None  # days after sowing to analyze minirhizotron data
    )
    MEASUREMENT_DAYS_MINIRHIZOTRON: list[str] = (
        None  # days after sowing when minirhizotron measurements were taken
    )

    # Growth period
    START_DATE: str = "19941106"
    END_DATE: str = "19950809"

    # File paths (string or Path; leave empty and set externally)
    PATH_TO_SOIL_TEMPERATURE_DATA: str | Path = ""
    PATH_TO_PLANT_PARAMETERS: str | Path = ""
    PATH_TO_SOIL_CORE_DATA: str | Path = ""
    PATH_TO_MINIRHIZOTRON_DATA_1: str | Path = ""
    PATH_TO_MINIRHIZOTRON_DATA_2: str | Path = ""
    PATH_TO_MINIRHIZOTRON_DATA_3: str | Path = ""
    OUTPUT_DIRECTORY: str | Path = DEFAULT_OUTPUT_DIR

    # Simulation parameters
    SIMULATION_TIME: int = 277  # days
    TIME_STEP: int = 1  # days
    N_SIMULATION_RUNS: int = 3  # number of simulation runs
    BASE_SEED_PLANT: int = 34567
    BASE_SEED: int = 56789
    BASE_SEED_SOIL_CORE: int = 12345


RUTHE_1994_95 = RutheConfig(
    START_DATE="19941106",
    END_DATE="19950809",
    SIMULATION_TIME=227,  # day of last measurement, total growth period actually 277 days
    PATH_TO_SOIL_TEMPERATURE_DATA=DATA_RUTHE_PROCESSED / "bodentemp_95.pkl",
    PATH_TO_PLANT_PARAMETERS=DATA_RUTHE_RAW / "wheat_1995_org1.xml",
    PATH_TO_SOIL_CORE_DATA=DATA_RUTHE_RAW / "SoilCores1995.csv",
    PATH_TO_MINIRHIZOTRON_DATA_1=DATA_RUTHE_RAW / "MiniRhizoDV260495.csv",
    PATH_TO_MINIRHIZOTRON_DATA_2=DATA_RUTHE_RAW / "MiniRhizoDV150595.csv",
    PATH_TO_MINIRHIZOTRON_DATA_3=DATA_RUTHE_RAW / "MiniRhizoDV130695.csv",
    DAYS_TO_ANALYZE_MINIRHIZOTRON=[172, 191, 220],
    MEASUREMENT_DAYS_MINIRHIZOTRON=[
        "26$^{\mathrm{th}}$ April, 1995",
        "15$^{\mathrm{th}}$ May, 1995",
        "13$^{\mathrm{th}}$ June, 1995",
    ],
    MEASUREMENT_DAY_SOIL_CORE="20$^{\mathrm{th}}$ June, 1995",
)


RUTHE_1995_96 = RutheConfig(
    START_DATE="19951102",
    END_DATE="19960827",
    SIMULATION_TIME=238,  # day of last measurement, total growth period actually 300 days
    PATH_TO_SOIL_TEMPERATURE_DATA=DATA_RUTHE_PROCESSED / "bodentemp_96.pkl",
    PATH_TO_PLANT_PARAMETERS=DATA_RUTHE_RAW / "wheat_1996_org1.xml",
    PATH_TO_SOIL_CORE_DATA=DATA_RUTHE_RAW / "SoilCores1996.csv",
    PATH_TO_MINIRHIZOTRON_DATA_1=DATA_RUTHE_RAW / "MiniRhizo1996.csv",
    DAYS_TO_ANALYZE_MINIRHIZOTRON=[181, 195, 209, 230],
    MEASUREMENT_DAYS_MINIRHIZOTRON=[
        "30$^{\mathrm{th}}$ April, 1996",
        "14$^{\mathrm{th}}$ May, 1996",
        "28$^{\mathrm{th}}$ May, 1996",
        "18$^{\mathrm{th}}$ June, 1996",
    ],
    MEASUREMENT_DAY_SOIL_CORE="26$^{\mathrm{th}}$ June, 1996",
)


RUTHE_1996_97 = RutheConfig(
    START_DATE="19961105",
    END_DATE="19970616",
    SIMULATION_TIME=224,
    PATH_TO_SOIL_TEMPERATURE_DATA=DATA_RUTHE_PROCESSED / "bodentemp_97.pkl",
    PATH_TO_PLANT_PARAMETERS=DATA_RUTHE_RAW / "wheat_1997_org1.xml",
    PATH_TO_SOIL_CORE_DATA=DATA_RUTHE_RAW / "SoilCores1997.csv",
    PATH_TO_MINIRHIZOTRON_DATA_1=DATA_RUTHE_RAW / "MiniRhizo1997.csv",
    DAYS_TO_ANALYZE_MINIRHIZOTRON=[162, 193, 207, 227, 240],
    MEASUREMENT_DAYS_MINIRHIZOTRON=[
        "15$^{\mathrm{th}}$ April, 1997",
        "16$^{\mathrm{th}}$ May, 1997",
        "30$^{\mathrm{th}}$ May, 1997",
        "19$^{\mathrm{th}}$ June, 1997",
        "02$^{\mathrm{th}}$ June, 1997",
    ],
    MEASUREMENT_DAY_SOIL_CORE="16$^{\mathrm{th}}$ June, 1997",
)
