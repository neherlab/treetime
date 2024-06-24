class AutoRepr(object):
  """
  Inherit from this to get JSON-like output when printing a class instance
  """

  def __repr__(self):
    items = ("%s: %r" % (k, v) for k, v in self.__dict__.items() if k != '__orig_class__')
    return f"{self.__class__.__name__} {{ {'; '.join(items)} }}"
