try:
    from enum import StrEnum
except ImportError:  # Python < 3.11
    from enum import Enum

    class StrEnum(str, Enum):
        """Minimal StrEnum backport for Python 3.9/3.10."""

        def __str__(self) -> str:
            return str(self.value)


class FvsVariant(StrEnum):
    """Enumeration of supported FVS Variants."""

    ALASKA = "AK"  # Southeast Alaska and Coastal British Columbia
    BLUE_MOUNTAINS = "BM"  # Blue Mountains
    INLAND_CALIFORNIA = "CA"  # Inland California and Southern Cascades (ICASCA)
    CENTRAL_IDAHO = "CI"  # Central Idaho
    CENTRAL_ROCKIES = "CR"  # Central Rockies
    CENTRAL_STATES = "CS"  # Central States
    EASTERN_CASCADES = "EC"  # East Cascades
    EASTERN_MONTANA = "EM"  # Eastern Montana
    INLAND_EMPIRE = "IE"  # Inland Empire
    KOOTENAI = "KT"  # Kootenai, Kaniksu, and Tally Lake (KooKanTL)
    LAKE_STATES = "LS"  # Lake States
    KLAMATH_MOUNTAINS = "NC"  # Klamath Mountains (and northern California)
    NORTHEAST_US = "NE"  # Northeastern US
    ORGANON_SOUTHWEST = "OC"  # Organon Southwest
    ORGANON_PACIFIC = "OP"  # Organon Pacific Northwest
    PACIFIC_COAST = "PN"  # Pacific Northwest Coast
    SOUTHERN_US = "SN"  # Southern US
    SOUTHERN_OREGON = (
        "SO"  # South Central Oregon and Northeast California (SORNEC)
    )
    TETONS = "TT"  # Tetons
    UTAH = "UT"  # Utah
    WESTERN_CASCADES = "WC"  # Westside Cascades
    WESTERN_SIERRAS = "WS"  # Western Sierra Nevada

    # fvs-modern additions (not in upstream USDA)
    ACADIAN = "ACD"  # Acadian (Maine, New Hampshire, Vermont)
    BRITISH_COLUMBIA = "BC"  # British Columbia
    ONTARIO = "ON"  # Ontario


class FvsAttributeAccessor(StrEnum):
    """How an FVS attribute is to be accessed."""

    GET = "get"
    SET = "set"
