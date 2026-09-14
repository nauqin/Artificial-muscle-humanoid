"""환경이 맞는지 먼저 보고, 아니면 어떻게 고치는지 알려준다.

**의존성 없이 표준 라이브러리만 쓴다.** numpy조차 import 하지 않는다 — 환경이
틀렸을 때 numpy가 없어 raw 트레이스백이 먼저 뜨면 안내가 묻히기 때문이다.

각 실행 스크립트 맨 위에서 `import _env` 만 하면 된다.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

MESSAGE = """
[환경 오류] 이 스크립트는 conda `{env}` 환경에서 돌아간다.

  지금 쓰는 파이썬 : {python}
  없는 것          : {missing}

고치는 방법은 둘 중 하나다.

  1) 어느 폴더에서든 한 방:
       {here}/run.sh {script} <인자…>

  2) 환경을 켜 두고 여러 명령을 이어서:
       conda activate {env}
       cd {here}
       python {script} <인자…>

※ VSCode의 ▶ 버튼이나 `/usr/local/bin/python3` 로는 안 된다 — 그쪽 파이썬에는
   OpenSim이 없다. 인터프리터를 `{env}` 로 바꾸거나 위 방법을 쓴다.
"""


def require(*modules, env="osim"):
    """필요한 모듈이 없으면 안내를 찍고 종료한다."""
    missing = []
    for name in modules:
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    if not missing:
        return
    print(MESSAGE.format(env=env, python=sys.executable, missing=", ".join(missing),
                         here=HERE, script=os.path.basename(sys.argv[0]) or "batch.py"),
          file=sys.stderr)
    raise SystemExit(2)


# import 만으로 확인되게 한다.
require("numpy", "opensim")
