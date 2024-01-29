from cath_alphaflow import settings
import logging
import pytest
import oracledb as Oracle

from click.testing import CliRunner


LOG = logging.getLogger(__name__)


class MockCursorIterator:
    """Iterator class"""

    def __init__(self, curs):
        self._index = 0
        self._curs = curs
        self._rows = curs.fetchmany()

    def __next__(self):
        if self._index < len(self._rows):
            result = self._rows[self._index]
            self._index += 1
            return result
        raise StopIteration


class MockCursor:
    def __iter__(self):
        """Returns the Iterator object"""
        return MockCursorIterator(self)

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def __init__(self, *args, **kwargs):
        self.fetched = False
        self.rows = None
        self.description = [("id",), ("name",)]

    def callfunc(self, *args):
        return self

    @staticmethod
    def var(*args):
        pass

    @staticmethod
    def execute(query, *args):
        pass

    @staticmethod
    def close():
        pass

    def fetchmany(self, *args):
        if not self.fetched:
            self.rows = self.fetch_rows()
            self.fetched = True

        return next(self.rows)

    @staticmethod
    def fetchall():
        return [
            (1, "TEST A"),
            (2, "TEST B"),
            (3, "TEST C"),
        ]

    def fetch_rows(self):
        yield self.fetchall()
        yield None


class MockConnect:
    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, type, value, traceback):
        pass

    @staticmethod
    def cursor():
        return MockCursor()

    @staticmethod
    def close():
        pass


@pytest.fixture
def mock_connection(monkeypatch, mock_settings):
    def mock_connect(*args, **kwargs):
        return MockConnect()

    monkeypatch.setattr(Oracle, "connect", mock_connect)


@pytest.fixture
def create_mock_query(monkeypatch, mock_connection):
    def _create_mock_query(description, rows):
        def get_mock_cursor(*args, **kwargs):
            mock_cursor = MockCursor()
            mock_cursor.description = description
            mock_cursor.fetchall = lambda: rows
            return mock_cursor

        monkeypatch.setattr(MockConnect, "cursor", get_mock_cursor)

    yield _create_mock_query

    # cleanup


@pytest.fixture
def mock_settings(monkeypatch):
    def mock_get_default_settings(*args, **kwargs):
        return settings.TestSettings()

    monkeypatch.setattr(settings, "get_default_settings", mock_get_default_settings)


@pytest.fixture
def create_cli_runner(monkeypatch):
    def _create_cli_runner(**kwargs):
        def mock_get_default_settings():
            _settings = settings.TestSettings()
            for key, val in kwargs.items():
                LOG.info(f"overriding local tests setting: {key}={val}")
                print(f"overriding local tests setting: {key}={val}")
                setattr(_settings, key, val)
            return _settings

        monkeypatch.setattr(settings, "get_default_settings", mock_get_default_settings)

        cli_runner = CliRunner()
        return cli_runner

    yield _create_cli_runner

    # cleanup
