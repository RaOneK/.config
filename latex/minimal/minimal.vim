set nocompatible
let &runtimepath  = '~/.local/share/nvim/site/pack/packer/start/vimtex,' . &runtimepath
let &runtimepath .= ',~/.local/share/nvim/site/pack/packer/start/vimtex/after'
filetype plugin indent on
syntax enable
" Add relevant options and VimTeX configuration below.
let g:vimtex_compiler_latexmk = {'continuous': 1}
