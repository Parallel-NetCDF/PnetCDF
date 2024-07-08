divert(`-1')
# list_len((item_1, item_2, ..., item_n))
#   parenthesized list, simple version
define(`list_len', `_list_len($@, 0)')`'dnl
define(`_list_len',`ifelse(`$1', `()', `$2', `$0((shift$1), incr(`$2'))')')`'dnl
divert`'dnl