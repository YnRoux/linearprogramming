# Matrix-based Simplex Method in MATLAB

처음 프로그래밍을 배우던 시기에 진행한 선형계획법 팀 프로젝트를 다시 정리한 저장소입니다. 팀 과제로 진행되었지만 문제 선정, 알고리즘 설계, MATLAB 솔버 구현 및 최종 발표는 Junha Ryu ([@YnRoux](https://github.com/YnRoux))가 단독으로 담당했습니다. 여러 실험 파일에 흩어져 있던 코드를 하나의 함수로 통합하고, 행 연산을 직접 반복하던 구현을 **2단계 revised simplex method**로 다시 작성했습니다.

이 저장소의 목적은 MATLAB 내장 함수의 대체품을 만드는 것이 아니라, 심플렉스법에서 기저행렬이 어떻게 사용되는지를 코드로 확인하는 것입니다.

## 작성자 및 기여

**Junha Ryu ([@YnRoux](https://github.com/YnRoux))**

- 프로젝트 주제와 문제 선정
- 선형계획법 솔버의 알고리즘 설계
- 기존 MATLAB 코드 전체 구현
- 팀 프로젝트 최종 발표
- revised simplex 기반 재설계 및 코드 정리
- 예제 통합, 자동 테스트와 문서 작성

## 지원하는 문제

기본 호출은 다음 문제를 풉니다.

```text
maximize    c' x
subject to  A x <= b
            x >= 0
```

옵션을 사용하면 다음 기능도 지원합니다.

- 최대화와 최소화
- `<=`, `>=`, `=` 제약의 혼합
- 음수 우변의 자동 정규화
- 비음수 제약이 없는 자유변수
- Phase I을 통한 초기 실행가능기저 탐색
- 실행 불가능 및 무한해 판정
- Bland 규칙을 이용한 피벗 동률 처리
- 중복 등식 제약 제거

일반적인 유한 하한·상한은 아직 직접 지원하지 않습니다. 필요한 경우 변수 치환 또는 추가 제약으로 표현해야 합니다.

## 사용법

`src` 폴더를 MATLAB 경로에 추가합니다.

```matlab
addpath('src')
```

예를 들어 다음 생산계획 문제를 풀 수 있습니다.

```matlab
c = [50; 40];
A = [2, 4;
     3, 1];
b = [80; 60];

[x, fval, exitflag, output] = simplexSolve(c, A, b)
```

결과는 다음과 같습니다.

```text
x = [16; 12]
fval = 1280
exitflag = 1
```

최소화, 혼합 제약 및 자유변수는 이름-값 옵션으로 지정합니다.

```matlab
[x, fval, exitflag, output] = simplexSolve(c, A, b, ...
    'ObjectiveSense', 'min', ...
    'ConstraintSense', [">="; "="; "<="], ...
    'FreeVariables', [false; true; false]);
```

`FreeVariables`에는 논리 벡터 대신 자유변수의 인덱스를 전달할 수도 있습니다.

## 행렬 기반 알고리즘

현재 기저변수의 열로 만든 행렬을 `B`라고 하면 기본해는 다음 연립방정식으로 계산합니다.

```text
x_B = B \ b
```

목적함수 계수 `c_B`로부터 쌍대가격을 구하고,

```text
y = B' \ c_B
```

각 변수의 감소비용을 계산합니다.

```text
r = c - A' y
```

최대화 문제에서 양의 감소비용을 가진 비기저변수가 없으면 현재 해가 최적입니다. 들어올 변수의 열을 `a_q`라고 하면 기저변수의 변화 방향은 다음과 같습니다.

```text
d = B \ a_q
```

`d_i > 0`인 행에 대해 `x_B(i)/d_i`의 최솟값을 구하여 나갈 변수를 선택합니다. 코드에서는 역행렬 `inv(B)`를 직접 만들지 않고 MATLAB의 역슬래시 연산자를 사용합니다.

초기 슬랙 기저를 바로 만들 수 없는 `>=` 및 등식 제약에는 인공변수를 추가합니다. Phase I에서 인공변수 합을 0으로 만들 수 있는지 검사한 다음, 인공변수를 제거하고 원래 목적함수로 Phase II를 실행합니다.

## 종료 상태

| `exitflag` | 의미 |
|---:|---|
| `1` | 최적해를 찾음 |
| `0` | 최대 반복 횟수에 도달 |
| `-2` | 실행 가능한 해가 없음 |
| `-3` | 목적함수가 무한히 개선될 수 있음 |
| `-4` | 특이 기저 등 수치 문제 발생 |

`output`에는 단계별 반복 횟수, 최종 기저, 표준형 해, 감소비용 및 제거된 중복 행이 들어 있습니다.

## 프로젝트 구조

```text
linearprogramming/
├─ README.md
├─ src/
│  └─ simplexSolve.m
├─ examples/
│  ├─ problemCatalog.m
│  └─ runExamples.m
└─ tests/
   └─ testSimplexSolve.m
```

- `src/simplexSolve.m`: 2단계 revised simplex 구현
- `examples/problemCatalog.m`: 기존 13개 예제를 하나의 데이터 카탈로그로 통합
- `examples/runExamples.m`: 모든 보존 예제를 실행
- `tests/testSimplexSolve.m`: 최적해, Phase I, 자유변수, 중복 제약, 무한해 및 실행 불가능 회귀 테스트

기존 `linpsol4`, `linpsol6`, `linpsol7`에서 같은 제약식 필드 이름을 재사용해 앞 제약이 사라지던 문제는 카탈로그에서 수정했습니다.

## 예제와 테스트 실행

전체 예제를 실행하려면 다음 명령을 사용합니다.

```matlab
cd examples
results = runExamples();
```

자동 테스트는 프로젝트 루트에서 실행합니다.

```matlab
results = runtests('tests');
table(results)
```

## 교육용 구현의 범위

이 구현은 작은 교육용 문제와 알고리즘 학습을 목표로 합니다. 대규모·희소·수치적으로 민감한 실제 최적화 문제에는 MATLAB Optimization Toolbox의 `linprog`처럼 전처리, 스케일링, 고급 피벗 전략을 갖춘 검증된 솔버를 사용하는 것이 적절합니다.
