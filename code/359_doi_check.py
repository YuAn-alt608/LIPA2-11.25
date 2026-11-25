#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
359 — 参考文献 DOI 程序化核对（Crossref）

对每条文献：有 DOI 的直接取 /works/{doi}；无 DOI 的按题名检索。
核对项：题名（词重叠）、期刊、卷、页、年份。输出 PASS / WARN 并打印 Crossref 原文，
不自动改写稿件——人工看过后由 360_r22_format.py 落盘。
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "ref-check/1.0 (mailto:yuan-zizheng@edu.ulisboa.pt)"}

REFS = [
    ("Adami 2025", "10.1016/j.xgen.2025.100979",
     "LINE-1 retrotransposons mediate cis-acting transcriptional control in human pluripotent "
     "stem cells and regulate early brain development", "Cell Genomics", "5", "100979"),
    ("Athanikar 2004", "10.1093/nar/gkh698",
     "A YY1-binding site is required for accurate human LINE-1 transcription initiation",
     "Nucleic Acids Res", "32", "3846-3855"),
    ("Bairoch 2018", "10.7171/jbt.18-2902-002",
     "The Cellosaurus, a cell-line knowledge resource", "J Biomol Tech", "29", "25-38"),
    ("Brown 1910", "10.1111/j.2044-8295.1910.tb00207.x",
     "Some experimental results in the correlation of mental abilities", None, "3", "296-322"),
    ("ENCODE 2012", "10.1038/nature11247",
     "An integrated encyclopedia of DNA elements in the human genome", "Nature", "489", "57-74"),
    ("Khan 2006", "10.1101/gr.4001406",
     "Molecular evolution and tempo of amplification of human LINE-1 retrotransposons since the "
     "origin of primates", "Genome Research", "16", "78-87"),
    ("Lanciano 2024", "10.1016/j.xgen.2024.100498",
     "Locus-level L1 DNA methylation profiling reveals the epigenetic and transcriptional "
     "interplay between L1s and their integration sites", "Cell Genomics", "4", "100498"),
    ("Liu 2018", "10.1038/nature25179",
     "Selective silencing of euchromatic L1s revealed by genome-wide screens for L1 regulators",
     "Nature", "553", "228-232"),
    ("Luo 2020", "10.1093/nar/gkz1062",
     "New developments on the Encyclopedia of DNA Elements (ENCODE) data portal",
     "Nucleic Acids Res", "48", "D882-D889"),
    ("Spearman 1904", "10.2307/1412159",
     "The proof and measurement of association between two things", "Am J Psychol", "15", "72-101"),
    ("Spearman 1910", "10.1111/j.2044-8295.1910.tb00206.x",
     "Correlation calculated from faulty data", None, "3", "271-295"),
    ("Stergachis 2020", "10.1126/science.aaz1646",
     "Single-molecule regulatory architectures captured by chromatin fiber sequencing",
     "Science", "368", "1449-1454"),
    ("Sun 2018", "10.1073/pnas.1722565115",
     "Transcription factor profiling reveals molecular choreography and key regulators of human "
     "retrotransposon expression", "Proc Natl Acad Sci USA", "115", "E5526-E5535"),
    ("Vollger 2026 (v3)", "10.1101/2024.06.14.599122",
     "Somatic epimutations cap genetic determinism in the human diploid chromatin epigenome",
     "bioRxiv", None, None),
]


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def norm(s):
    s = (s or "").lower()
    s = s.replace("–", "-").replace("—", "-")
    return re.sub(r"[^a-z0-9 ]+", " ", s)


def overlap(a, b):
    A, B = set(norm(a).split()), set(norm(b).split())
    return len(A & B) / max(1, len(A))


def show(tag, m, doi):
    t = (m.get("title") or [""])[0]
    print(f"    CR title : {t}")
    print(f"    CR DOI   : {m.get('DOI')}")
    print(f"    CR cont  : {(m.get('container-title') or [''])[0][:60]} | vol {m.get('volume')}"
          f" | page {m.get('page')} | year {m.get('issued', {}).get('date-parts', [['?']])[0][0]}")


def main():
    fails = []
    for tag, doi, title, journal, vol, page in REFS:
        print(f"\n=== {tag}")
        try:
            if doi:
                m = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
            else:
                q = urllib.parse.urlencode({"query.bibliographic": title, "rows": 5})
                got = get("https://api.crossref.org/works?" + q)["message"]["items"]
                if not got:
                    print("    !! no crossref hit")
                    fails.append((tag, "no-hit"))
                    continue
                m = got[0]
                # 取标题重合度最高的一条
                m = max(got, key=lambda x: overlap(title, (x.get("title") or [""])[0]))
            show(tag, m, doi)
            ov = overlap(title, (m.get("title") or [""])[0])
            ok_title = ov >= 0.75
            ok_journal = True if journal is None else \
                journal.lower()[:12] in ((m.get("container-title") or [""])[0] or "").lower()
            ok_vol = True if vol is None else str(m.get("volume")) == str(vol)
            ok_page = True if page is None else \
                (m.get("page") or "").replace("–", "-").replace(" ", "") == page.replace("-", "-")
            status = "PASS" if all([ok_title, ok_journal, ok_vol, ok_page]) else "WARN"
            print(f"    -> {status}  title {ov:.2f}  journal {ok_journal}  vol {ok_vol}  page {ok_page}")
            if status != "PASS":
                fails.append((tag, f"title={ov:.2f} j={ok_journal} v={ok_vol} p={ok_page}"))
            if not doi:
                print(f"    FOUND DOI: {m.get('DOI')}")
        except Exception as e:
            print(f"    !! {e}")
            fails.append((tag, str(e)))
        time.sleep(0.4)

    print("\n---- summary ----")
    for t, r in fails:
        print(f"  {t}: {r}")
    print("fails:", len(fails))


if __name__ == "__main__":
    main()
