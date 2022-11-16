class BaseError(BaseException):
    pass


class ParseError(BaseError):
    pass


class UsageError(BaseError):
    pass


class MultipleModelsError(BaseError):
    pass


class MultipleChainsError(BaseError):
    pass


class ChoppingError(BaseError):
    pass


class NoMatchingResiduesError(BaseError):
    pass
