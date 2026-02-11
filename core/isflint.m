function isf = isflint(m)
isf = (abs(m) <= flintmax && m == floor(m));