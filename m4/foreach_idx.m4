divert(`-1')
# foreach_idx(x, idx, (item_1, item_2, ..., item_n), stmt)
#   parenthesized list, simple version
define(`foreach_idx', `pushdef(`$1')pushdef(`$2')_foreach_idx($@,0)popdef(`$2')popdef(`$1')')
define(`_arg1', `$1')
define(`_foreach_idx', `ifelse(`$3', `()', `',`define(`$1', _arg1$3)define(`$2', `$5')$4`'$0(`$1', `$2', (shift$3), `$4',incr($5))')')
divert`'dnl