class BaseError(BaseException):
    pass


class ParseError(BaseError):
    pass


class ArgumentError(BaseError):
    pass


class NoMatchingFragmentError(BaseError):
    pass


class MultipleModelsError(BaseError):
    pass


class MultipleChainsError(BaseError):
    pass


class ChoppingError(BaseError):
    pass


class NoMatchingResiduesError(BaseError):
    pass
