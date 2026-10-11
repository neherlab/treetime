import type { PortMessage as AddonPortMessage, PortReply as AddonPortReply } from "@neherlab/app-napi";
import { type PortMessage, type PortReply, portReplySchema } from "@neherlab/app-ui/port";

describe("port_types of the addon and the shared transport", () => {
  test("every reply of the addon is a reply of the shared transport, and every message of the transport reaches the addon", () => {
    expectTypeOf<AddonPortReply>().toExtend<PortReply>();
    expectTypeOf<PortMessage>().toExtend<AddonPortMessage>();
  });

  test("the replies of the addon pass the reply schema of the shared transport", () => {
    const replies: AddonPortReply[] = [
      { kind: "head", seq: 0, status: 200, headers: [{ name: "content-type", value: "text/plain" }] },
      { kind: "chunk", seq: 0, data: new Uint8Array([104, 105]) },
      { kind: "end", seq: 0 },
      { kind: "reset", seq: 1, message: "When reading the response body: gone" },
    ];

    expect(replies.map((reply) => portReplySchema.parse(reply))).toStrictEqual(replies);
  });
});
