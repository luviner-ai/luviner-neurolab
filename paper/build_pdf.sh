#!/bin/bash
# Rebuild preprint.pdf from preprint.md.
# Requires pandoc, plus one HTML-to-PDF renderer: headless Chrome or
# Chromium (preferred, it renders the CSS as designed) or wkhtmltopdf.
set -euo pipefail
cd "$(dirname "$0")"

command -v pandoc >/dev/null || { echo "pandoc not found: install it and retry" >&2; exit 1; }

mkdir -p build
cat > build/style.css << 'CSS'
body{font-family:Georgia,"Times New Roman",serif;font-size:11pt;line-height:1.45;max-width:17cm;margin:2cm auto;color:#111}
h1{font-size:18pt;line-height:1.25} h2{font-size:14pt;margin-top:1.6em;border-bottom:1px solid #999}
h3{font-size:12pt} h4{font-size:11pt;font-style:italic}
table{border-collapse:collapse;font-size:9pt;margin:1em 0} th,td{border:1px solid #999;padding:3px 6px;vertical-align:top}
img{max-width:100%;display:block;margin:1em auto} code{font-size:9.5pt}
blockquote{border-left:3px solid #c00;padding-left:.8em;color:#900}
@page{size:A4;margin:2cm}
CSS
cp preprint.md build/preprint.md
pandoc build/preprint.md --from gfm --to html5 --standalone --embed-resources --css style.css \
  --metadata title=" " --resource-path=.:figures -o build/preprint.html

CHROME=""
for c in "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
         "/Applications/Chromium.app/Contents/MacOS/Chromium" \
         "$(command -v google-chrome || true)" \
         "$(command -v google-chrome-stable || true)" \
         "$(command -v chromium || true)" \
         "$(command -v chromium-browser || true)"; do
  [ -n "$c" ] && [ -x "$c" ] && { CHROME="$c"; break; }
done

if [ -n "$CHROME" ]; then
  "$CHROME" --headless=new --disable-gpu --no-pdf-header-footer \
    --print-to-pdf="$PWD/preprint.pdf" "file://$PWD/build/preprint.html" >/dev/null 2>&1
elif command -v wkhtmltopdf >/dev/null; then
  echo "note: falling back to wkhtmltopdf; page breaks may differ from the shipped PDF" >&2
  wkhtmltopdf --enable-local-file-access -q build/preprint.html preprint.pdf
else
  cat >&2 <<'MSG'
No HTML-to-PDF renderer found. Install ONE of:
  - Google Chrome or Chromium (preferred; renders the CSS as designed)
  - wkhtmltopdf
build/preprint.html has been written and can be printed to PDF by hand.
MSG
  exit 1
fi
rm -rf build
echo "wrote preprint.pdf ($(du -h preprint.pdf | cut -f1))"
