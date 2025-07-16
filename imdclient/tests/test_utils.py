# tests/test_utils.py
import pytest
from imdclient.utils import parse_host_port


@pytest.mark.parametrize(
    "server_address,host_port",
    [
        ("imd://localhost:8888", ("localhost", 8888)),
        ("imd://example.com:12345", ("example.com", 12345)),
    ],
)
def test_parse_host_port_valid(server_address, host_port):
    assert parse_host_port(server_address) == host_port


@pytest.mark.parametrize(
    "server_address",
    [
        "",  # empty
        "http://localhost:80",  # wrong protocol prefix
        "imd://",  # missing host and port
        "imd://localhost:",  # missing port
        "imd://:8080",  # missing host
        "imd://host:notaport",  # port not integer
        "imd://host:80:90",  # too many segments
    ],
)
def test_parse_host_port_invalid(server_address):
    with pytest.raises(ValueError):
        parse_host_port(server_address)
